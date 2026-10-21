// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/ProteoformDatabase.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <set>
#include <unordered_map>

namespace OpenMS
{
namespace
{
// the residues a decoy permutation may move (any other, e.g. X or U, stays in place)
const std::string movable_residues = "ACDEFGHIKLMNPQRSTVWY";
// the half cystine of UniPEFF's disulfide bonds: used only through the \DisulfideBond pairs
const String half_cystine = "MOD:00798";
// at most this many modifications of one residue form its unit (2^n - 1 choices)
constexpr Size max_mods_per_residue = 8;

// PEFF:00010xx accessions of \Processed (UniPEFF): the kept chain and the removed segments
const String processed_mature = "PEFF:0001020";
const std::set<String> processed_removed {"PEFF:0001035", "PEFF:0001021", "PEFF:0001022", "PEFF:0001034"};

// a key of a PEFF description: an unescaped backslash followed by [A-Za-z0-9_]+ and '='
bool isKeyChar(char c)
{
  return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_';
}

bool isEscapable(char c)
{
  return c == '\\' || c == '|' || c == '(' || c == ')';
}

// the \Key=value pairs of a PEFF description (values trimmed, still escaped); PEFFFile.cpp:2049-2078 upstream
std::map<String, String> parseKeys(const String& description)
{
  std::map<String, String> keys;
  String key, value;
  bool in_value = false;
  auto flush = [&]() {
    if (! in_value) return;
    value.trim();
    if (keys.find(key) == keys.end()) keys[key] = value;
    value.clear();
  };
  for (Size i = 0; i < description.size(); i++)
  {
    const char c = description[i];
    if (c == '\\' && i + 1 < description.size())
    {
      if (isEscapable(description[i + 1]))
      {
        if (in_value) value += description.substr(i, 2);
        i++;
        continue;
      }
      Size j = i + 1;
      while (j < description.size() && isKeyChar(description[j])) j++;
      if (j > i + 1 && j < description.size() && description[j] == '=')
      {
        flush();
        key = description.substr(i + 1, j - i - 1);
        in_value = true;
        i = j;
        continue;
      }
    }
    if (in_value) value += c;
  }
  flush();
  return keys;
}

// the top-level parenthesised tuples of a value, by depth, skipping escaped characters (PEFFFile.cpp:2699-2767 upstream); a value
// without parentheses is one tuple
std::vector<String> splitTuples(const String& value)
{
  std::vector<String> tuples;
  if (value.empty() || value[0] != '(')
  {
    if (! value.empty()) tuples.push_back(value);
    return tuples;
  }
  int depth = 0;
  String current;
  for (Size i = 0; i < value.size(); i++)
  {
    const char c = value[i];
    if (c == '\\' && i + 1 < value.size() && isEscapable(value[i + 1]))
    {
      if (depth > 0) current += value.substr(i, 2);
      i++;
      continue;
    }
    if (c == '(')
    {
      if (depth > 0) current += c;
      depth++;
      continue;
    }
    if (c == ')')
    {
      depth--;
      if (depth == 0)
      {
        tuples.push_back(current);
        current.clear();
        continue;
      }
      if (depth < 0) depth = 0;
    }
    if (depth > 0) current += c;
  }
  return tuples;
}

String unescape(const String& s)
{
  String out;
  out.reserve(s.size());
  for (Size i = 0; i < s.size(); i++)
  {
    if (s[i] == '\\' && i + 1 < s.size() && isEscapable(s[i + 1]))
    {
      out += s[i + 1];
      i++;
      continue;
    }
    out += s[i];
  }
  return out;
}

// the fields of a tuple, split at unescaped '|' and unescaped (PEFFFile.cpp:35-112 upstream)
std::vector<String> splitFields(const String& tuple)
{
  std::vector<String> fields;
  String current;
  for (Size i = 0; i < tuple.size(); i++)
  {
    if (tuple[i] == '\\' && i + 1 < tuple.size() && isEscapable(tuple[i + 1]))
    {
      current += tuple.substr(i, 2);
      i++;
      continue;
    }
    if (tuple[i] == '|')
    {
      fields.push_back(unescape(current));
      current.clear();
      continue;
    }
    current += tuple[i];
  }
  fields.push_back(unescape(current));
  return fields;
}

// the Option B identifier "id:" of the first field: stripped and returned (-1: none)
long stripIdentifier(String& field)
{
  const Size colon = field.find(':');
  if (colon == String::npos || colon == 0) return -1;
  for (Size i = 0; i < colon; i++)
  {
    if (field[i] < '0' || field[i] > '9') return -1;
  }
  const long id = std::stol(field.substr(0, colon));
  field = field.substr(colon + 1);
  return id;
}

bool toInt(const String& s, int& value)
{
  String t = s;
  t.trim();
  if (t.empty()) return false;
  for (char c : t)
  {
    if (c < '0' || c > '9') return false;
  }
  value = std::stoi(t);
  return true;
}

// the accession of an identifier: ACC of sp:ACC, sp|ACC|ID or ACC|..., else the identifier
String accessionOf(const String& identifier)
{
  std::vector<String> parts;
  identifier.split('|', parts);
  if (parts.size() >= 3) return parts[1];
  const Size colon = identifier.find(':');
  if (colon != String::npos) return identifier.substr(colon + 1);
  return identifier;
}

// a modification of a residue of a protein (1-based positions)
struct SiteModification
{
  int position = 0;
  double mass = 0;
  String accession, name;
  int sites = 1;  // the sites of the residue it occupies (FixedModification::sites)
  int bond = -1;  // a half cystine: the other one of its bond
};

// The TermSpec and Origin of the PSI-MOD terms (the bundled PSI-MOD.obo), read once. ModificationsDB keeps neither for a term it maps to
// a Unimod record: MOD:00069 (N-palmitoyl-L-cysteine, N-term) is "Palmitoyl (C)" there, anywhere.
struct PsiModTerm
{
  String term_spec, origin;
};

const std::map<String, PsiModTerm>& psiModTerms()
{
  static const std::map<String, PsiModTerm> terms = [] {
    std::map<String, PsiModTerm> t;
    String path;
    try
    {
      path = File::find("CHEMISTRY/PSI-MOD.obo");
    }
    catch (const Exception::BaseException&)
    {
      return t;
    }
    std::ifstream in(path);
    std::string line;
    String id;
    // property_value: TermSpec: "N-term" xsd:string
    auto value = [&line](const std::string& key, String& out) {
      const std::string prefix = "property_value: " + key + ": \"";
      if (line.compare(0, prefix.size(), prefix) != 0) return;
      const Size end = line.find('"', prefix.size());
      if (end != std::string::npos) out = line.substr(prefix.size(), end - prefix.size());
    };
    while (std::getline(in, line))
    {
      if (! line.empty() && line[0] == '[')
      {
        id.clear();
        continue;
      }
      if (line.compare(0, 4, "id: ") == 0)
      {
        id = line.substr(4);
        id.trim();
        continue;
      }
      if (id.empty() || line.compare(0, 16, "property_value: ") != 0) continue;
      value("TermSpec", t[id].term_spec);
      value("Origin", t[id].origin);
    }
    return t;
  }();
  return terms;
}

// The sites of its residue a modification occupies (bits: 1 side chain, 2 alpha-amino, 4 alpha-carboxyl group), by PSI-MOD's TermSpec.
// Without one (a term the bundled PSI-MOD lacks, a Unimod accession): the side chain, and on a residue where a form starts (kind 2) or
// ends (kind 3) also that terminus.
int sitesOf(const String& accession, int kind)
{
  const auto& terms = psiModTerms();
  if (const auto it = terms.find(accession); it != terms.end())
  {
    if (it->second.term_spec == "N-term") return 2;
    if (it->second.term_spec == "C-term") return 4;
    if (it->second.term_spec == "none") return 1;
  }
  return kind == 2 ? 3 : (kind == 3 ? 5 : 1);
}

// whether a modification can sit on a residue: its origin (PSI-MOD Origin; the residues of the Unimod records) is the residue or any
// residue (X); unknown: yes
bool residueFits(const String& accession, char residue)
{
  if (accession.hasPrefix("MOD:"))
  {
    const auto& terms = psiModTerms();
    const auto it = terms.find(accession);
    if (it == terms.end() || it->second.origin.empty()) return true;
    std::vector<String> origins;
    it->second.origin.split(',', origins);
    for (auto o : origins)
    {
      o.trim();
      if (o == "X" || (o.size() == 1 && o[0] == residue)) return true;
    }
    return false;
  }
  if (accession.hasPrefix("UNIMOD:"))
  {
    const String name = "UniMod:" + accession.substr(7);
    const auto* mdb = ModificationsDB::getInstance();
    if (! mdb->has(name)) return true;
    std::set<const ResidueModification*> mods;
    mdb->searchModifications(mods, name);
    if (mods.empty()) return true;
    for (const auto* m : mods)
    {
      if (m->getOrigin() == 'X' || m->getOrigin() == residue) return true;
    }
    return false;
  }
  return true;
}

// the identifier without its decoy prefix (DECOY_ or DECOY)
String withoutDecoyPrefix(const String& identifier)
{
  if (identifier.hasPrefix("DECOY_")) return identifier.substr(6);
  if (identifier.hasPrefix("DECOY")) return identifier.substr(5);
  return identifier;
}
} // namespace

std::uint64_t ProteoformDatabase::fnv1a64(const String& s)
{
  std::uint64_t h = 14695981039346656037ULL;
  for (const char c : s)
  {
    h ^= (std::uint64_t)(unsigned char)c;
    h *= 1099511628211ULL;
  }
  return h;
}

String ProteoformDatabase::proformaTag(const String& accession, double mass)
{
  return ProForma::formatMass(mass) + "|" + accession;
}

bool ProteoformDatabase::isPEFF(const String& file)
{
  if (String(file).toLower().hasSuffix(".peff")) return true;
  std::ifstream in(file);
  String line;
  if (! std::getline(in, line)) return false;
  return line.hasPrefix("# PEFF");
}

ProteoformDatabase ProteoformDatabase::fromFASTA(const std::vector<FASTAFile::FASTAEntry>& fasta_entries)
{
  ProteoformDatabase db;
  db.proteins.resize(fasta_entries.size());
  db.entries.resize(fasta_entries.size());
  // the partner of a protein by name, as FLASHTnT found it: DECOY_ or DECOY + the name of a target, the first entry of a name
  std::unordered_map<std::string, Size> index_of;
  index_of.reserve(fasta_entries.size());
  for (Size i = 0; i < fasta_entries.size(); i++) index_of.emplace(fasta_entries[i].identifier, i);
  for (Size i = 0; i < fasta_entries.size(); i++)
  {
    const auto& fe = fasta_entries[i];
    auto& p = db.proteins[i];
    p.identifier = fe.identifier;
    p.description = fe.description;
    p.sequence = fe.sequence;
    p.decoy = fe.identifier.hasPrefix("DECOY");
    p.entries = {i};
    p.forms = {{1, (int)fe.sequence.size()}};
    std::vector<std::string> names;
    if (fe.identifier.hasPrefix("DECOY_")) { names.push_back(fe.identifier.substr(6)); }
    else if (fe.identifier.hasPrefix("DECOY")) { names.push_back(fe.identifier.substr(5)); }
    else { names = {"DECOY_" + fe.identifier, "DECOY" + fe.identifier}; }
    for (const auto& name : names)
    {
      const auto it = index_of.find(name);
      if (it == index_of.end()) continue;
      p.partner = (long)it->second;
      break;
    }
    auto& e = db.entries[i];
    e.protein = i;
    e.offset = 0;
    e.length = (int)fe.sequence.size();
    e.id = fe.identifier;
  }
  db.finish_();
  return db;
}

void ProteoformDatabase::finish_()
{
  const Size n = entries.size();
  fixed_masses_.assign(n, {});
  occupancy_.assign(n, {});
  protein_of_.assign(n, 0);
  with_fixed_ = false;
  identity_ = entries.size() == proteins.size();
  for (Size e = 0; e < n; e++)
  {
    const auto& entry = entries[e];
    protein_of_[e] = entry.protein;
    identity_ = identity_ && entry.protein == e && entry.offset == 0 && entry.fixed.empty()
                && entry.length == (int)proteins[entry.protein].sequence.size();
    for (const auto& fm : entry.fixed)
    {
      with_fixed_ = true;
      auto& masses = fixed_masses_[e];
      auto& occupied = occupancy_[e];
      // the terminus of a kind-2 (3) modification, or the side chain of a kind-1 one, and the sites it occupies
      const std::uint8_t bit = (std::uint8_t)((fm.kind == 2 ? 2 : (fm.kind == 3 ? 4 : 1)) | (fm.sites & 7));
      if (! masses.empty() && masses.back().first == fm.position)
      {
        masses.back().second += fm.mass;
        occupied.back().second |= bit;
      }
      else
      {
        masses.emplace_back(fm.position, fm.mass);
        occupied.emplace_back(fm.position, bit);
      }
    }
  }
}

void ProteoformDatabase::loadPEFF(const String& file, int max_mods, const std::vector<String>& mod_masses,
                                  std::vector<FASTAFile::FASTAEntry>& fasta_entries)
{
  std::vector<FASTAFile::FASTAEntry> input;
  FASTAFile().load(file, input); // skips the '#' header lines of a PEFF file
  build(input, max_mods, mod_masses, fasta_entries);
}

void ProteoformDatabase::build(const std::vector<FASTAFile::FASTAEntry>& input, int max_mods, const std::vector<String>& mod_masses,
                               std::vector<FASTAFile::FASTAEntry>& fasta_entries)
{
  entries.clear();
  proteins.clear();
  fasta_entries.clear();
  statistics_ = LoadStatistics();
  max_mods = std::max(0, max_mods);

  // masses by accession given by the user: a mass, or "skip"
  std::map<String, double> given_masses;
  std::set<String> skipped;
  for (const auto& item : mod_masses)
  {
    String s = item;
    s.trim();
    if (s.empty()) continue;
    const Size eq = s.find('=');
    if (eq == String::npos)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "peff_mod_masses: expected accession=mass or accession=skip", s);
    }
    String acc = s.substr(0, eq), value = s.substr(eq + 1);
    acc.trim();
    value.trim();
    if (String(value).toLower() == "skip") { skipped.insert(acc); }
    else { given_masses[acc] = value.toDouble(); }
  }
  // the mass of an accession on a residue (cached per accession and residue); false: none (or skipped)
  std::map<std::pair<String, char>, std::pair<bool, double>> mass_cache;
  auto massOf = [&](const String& acc, char residue, double& mass) {
    if (skipped.count(acc) > 0) return false;
    if (const auto it = given_masses.find(acc); it != given_masses.end())
    {
      mass = it->second;
      return true;
    }
    const auto key = std::make_pair(acc, residue);
    auto it = mass_cache.find(key);
    if (it == mass_cache.end())
    {
      std::pair<bool, double> result {false, 0.0};
      try
      {
        result.second = ProForma::getTagMass(acc, residue);
        result.first = std::abs(result.second) >= 1e-6;
      }
      catch (const Exception::BaseException&)
      {
      }
      it = mass_cache.emplace(key, result).first;
    }
    mass = it->second.second;
    return it->second.first;
  };
  const auto* rdb = ResidueDB::getInstance();
  auto residueMass = [&](char c) {
    if (movable_residues.find(c) == std::string::npos && c != 'U' && c != 'O') return 0.0;
    return rdb->getResidue(c)->getMonoWeight(Residue::Internal);
  };

  // A parsed protein: forms, units and the entries made of them. The entries of a decoy generated for it mirror them.
  struct Form
  {
    int start, end; // 1-based, inclusive
  };
  struct Unit
  {
    int position;                                        // its first residue (1-based), the order of units
    std::vector<std::vector<SiteModification>> choices;  // each a set of modifications (1-based positions)
  };
  struct ParsedProtein
  {
    String identifier, description, sequence, accession;
    bool decoy = false;
    std::vector<Form> forms;
    // per form, its entries: the fixed modifications of each (1-based positions, kind)
    std::vector<std::vector<std::vector<std::pair<SiteModification, int>>>> form_entries;
  };

  auto parseProtein = [&](const FASTAFile::FASTAEntry& fe, bool is_decoy) {
    ParsedProtein pp;
    pp.sequence = fe.sequence;
    pp.decoy = is_decoy;
    const int L = (int)fe.sequence.size();
    const auto keys = parseKeys(fe.description);
    const bool peff = keys.count("PName") > 0 || keys.count("ModResPsi") > 0 || keys.count("Processed") > 0 || keys.count("ID") > 0
                      || keys.count("ModResUnimod") > 0 || keys.count("NcbiTaxId") > 0;
    // the identifier sp|ACC|ID (as UniProt FASTA files) and a FASTA-like description; a decoy keeps its prefix
    String prefix, id = fe.identifier;
    if (is_decoy)
    {
      prefix = id.hasPrefix("DECOY_") ? "DECOY_" : "DECOY";
      id = id.substr(prefix.size());
    }
    if (peff)
    {
      const auto it = keys.find("ID");
      const Size colon = id.find(':');
      if (it != keys.end() && colon != String::npos) id = id.substr(0, colon) + "|" + id.substr(colon + 1) + "|" + unescape(it->second);
      String description;
      if (const auto pn = keys.find("PName"); pn != keys.end())
      {
        const auto names = splitTuples(pn->second);
        if (! names.empty()) description = unescape(names[0]);
      }
      auto append = [&](const String& key, const String& label) {
        const auto k = keys.find(key);
        if (k == keys.end()) return;
        String v = unescape(k->second);
        if (v.hasPrefix("(") && v.hasSuffix(")")) v = v.substr(1, v.size() - 2);
        description += (description.empty() ? "" : " ") + label + "=" + v;
      };
      append("TaxName", "OS");
      append("NcbiTaxId", "OX");
      append("GName", "GN");
      append("PE", "PE");
      append("SV", "SV");
      pp.description = description;
      if (const auto len = keys.find("Length"); len != keys.end())
      {
        int length = 0;
        if (toInt(len->second, length) && length != L) statistics_.length_mismatches++;
      }
    }
    else { pp.description = fe.description; }
    pp.identifier = prefix + id;
    pp.accession = accessionOf(id);

    // processed forms: the full length, every mature chain, the full length without a removed segment at its start (residue 1 or 2) or
    // at its end
    std::set<std::pair<int, int>> form_set {{1, L}};
    if (const auto pr = keys.find("Processed"); pr != keys.end())
    {
      for (const auto& tuple : splitTuples(pr->second))
      {
        auto fields = splitFields(tuple);
        if (fields.size() < 3) continue;
        stripIdentifier(fields[0]);
        int b = 0, e = 0;
        if (! toInt(fields[0], b) || ! toInt(fields[1], e) || b < 1 || e > L || b > e) continue;
        String acc = fields[2];
        acc.trim();
        if (acc == processed_mature) form_set.insert({b, e});
        else if (processed_removed.count(acc) > 0)
        {
          if (b <= 2 && e < L) form_set.insert({e + 1, L});
          else if (e == L && b > 1) form_set.insert({1, b - 1});
        }
      }
    }
    for (const auto& [b, e] : form_set) pp.forms.push_back({b, e});
    std::set<int> starts, ends;
    for (const auto& f : pp.forms)
    {
      starts.insert(f.start);
      ends.insert(f.end);
    }

    // modifications by residue; half cystines by their tuple identifier for the bonds
    std::map<int, std::vector<SiteModification>> site_mods;
    std::map<long, int> half_cystine_of; // tuple identifier -> position
    Size half_cystines = 0;
    for (const String key : {"ModResPsi", "ModResUnimod"})
    {
      const auto it = keys.find(key);
      if (it == keys.end()) continue;
      for (const auto& tuple : splitTuples(it->second))
      {
        auto fields = splitFields(tuple);
        if (fields.size() < 2) continue;
        const long tuple_id = stripIdentifier(fields[0]);
        int pos = 0;
        if (! toInt(fields[0], pos) || pos < 1 || pos > L) continue; // '?' and positions outside the sequence
        String acc = fields[1];
        acc.trim();
        if (acc.empty()) continue;
        const char residue = fe.sequence[pos - 1];
        if (! residueFits(acc, residue))
        {
          statistics_.wrong_residue[acc]++;
          continue;
        }
        if (acc == half_cystine)
        {
          half_cystines++;
          if (tuple_id >= 0) half_cystine_of[tuple_id] = pos;
          continue;
        }
        double mass = 0;
        if (! massOf(acc, residue, mass) || residueMass(residue) + mass <= 1.0)
        {
          if (! is_decoy) statistics_.unresolved[acc]++;
          continue;
        }
        auto& mods = site_mods[pos];
        if (mods.size() >= max_mods_per_residue) continue;
        mods.push_back({pos, mass, acc, fields.size() > 2 ? fields[2] : String()});
        mods.back().sites = sitesOf(acc, starts.count(pos) > 0 ? 2 : (ends.count(pos) > 0 ? 3 : 1));
        if (! is_decoy) statistics_.modification_sites++;
      }
    }
    // disulfide bonds: the pairs of half cystines of \DisulfideBond (UniPEFF -annotation_identifiers)
    std::vector<std::pair<int, int>> bonds;
    Size paired = 0;
    if (const auto it = keys.find("DisulfideBond"); it != keys.end())
    {
      for (const auto& tuple : splitTuples(it->second))
      {
        auto fields = splitFields(tuple);
        stripIdentifier(fields[0]);
        std::vector<String> ids;
        fields[0].split(',', ids);
        if (ids.size() != 2) continue;
        int a = 0, b = 0;
        if (! toInt(ids[0], a) || ! toInt(ids[1], b)) continue;
        const auto ia = half_cystine_of.find(a), ib = half_cystine_of.find(b);
        if (ia == half_cystine_of.end() || ib == half_cystine_of.end() || ia->second == ib->second) continue;
        bonds.emplace_back(std::min(ia->second, ib->second), std::max(ia->second, ib->second));
        paired += 2;
      }
    }
    double half_cystine_mass = 0;
    const bool with_bonds = ! bonds.empty() && massOf(half_cystine, 'C', half_cystine_mass);
    if (! is_decoy)
    {
      statistics_.unpaired_half_cystines += half_cystines - std::min(half_cystines, paired);
      if (with_bonds) statistics_.disulfide_bonds += bonds.size();
      statistics_.forms += pp.forms.size();
    }

    // the entries of each form: the form times every set of at most max_mods units, by size, then by unit (position); a unit is a residue
    // with a nonempty subset of its modifications (by bit mask) on distinct sites, or all disulfide bonds inside the form
    for (const auto& f : pp.forms)
    {
      std::vector<Unit> units;
      std::vector<int> unit_kind; // per unit: the kind of its modifications
      for (const auto& [pos, mods] : site_mods)
      {
        if (pos < f.start || pos > f.end) continue;
        // the N-terminal rule: a modification of a residue where some form starts is used only in forms that start there (kind 2); the
        // mirror for a residue where some form ends (kind 3)
        int kind = 1;
        if (starts.count(pos) > 0) kind = 2;
        else if (ends.count(pos) > 0) kind = 3;
        if (kind == 2 && pos != f.start) continue;
        if (kind == 3 && pos != f.end) continue;
        Unit u;
        u.position = pos;
        const Size n = mods.size();
        // the modifications of a choice sit on distinct sites of the residue: alternatives on one site (N6-acetyl or N6-succinyl) are never
        // combined, the lipoprotein pair (N-palmitoyl on the amino group, S-diacylglycerol on the side chain) is
        for (Size mask = 1; mask < ((Size)1 << n); mask++)
        {
          std::vector<SiteModification> choice;
          int sites = 0;
          for (Size m = 0; m < n && sites >= 0; m++)
          {
            if ((mask & ((Size)1 << m)) == 0) continue;
            if ((sites & mods[m].sites) != 0) sites = -1;
            else
            {
              sites |= mods[m].sites;
              choice.push_back(mods[m]);
            }
          }
          if (sites > 0) u.choices.push_back(std::move(choice));
        }
        units.push_back(std::move(u));
        unit_kind.push_back(kind);
      }
      if (with_bonds)
      {
        std::vector<SiteModification> all;
        for (const auto& [a, b] : bonds)
        {
          if (a < f.start || b > f.end) continue;
          all.push_back({a, half_cystine_mass, half_cystine, "half cystine", 1, b});
          all.push_back({b, half_cystine_mass, half_cystine, "half cystine", 1, a});
        }
        if (! all.empty())
        {
          Unit u;
          u.position = all.front().position;
          u.choices.push_back(std::move(all));
          units.push_back(std::move(u));
          unit_kind.push_back(1);
        }
      }
      // units by position (a disulfide unit at its first cysteine)
      std::vector<Size> order(units.size());
      for (Size i = 0; i < order.size(); i++) order[i] = i;
      std::stable_sort(order.begin(), order.end(), [&](Size a, Size b) { return units[a].position < units[b].position; });

      std::vector<std::vector<std::pair<SiteModification, int>>> form_entries;
      const int k_max = std::min(max_mods, (int)units.size());
      std::vector<Size> combination;
      for (int k = 0; k <= k_max; k++)
      {
        // combinations of k units in lexicographic order, each with every product of its units' choices
        combination.resize(k);
        for (int i = 0; i < k; i++) combination[i] = i;
        while (true)
        {
          std::vector<Size> choice(k, 0);
          while (true)
          {
            std::vector<std::pair<SiteModification, int>> fixed;
            for (int i = 0; i < k; i++)
            {
              const Size u = order[combination[i]];
              for (const auto& m : units[u].choices[choice[i]]) fixed.emplace_back(m, unit_kind[u]);
            }
            std::stable_sort(fixed.begin(), fixed.end(), [](const auto& a, const auto& b) { return a.first.position < b.first.position; });
            form_entries.push_back(std::move(fixed));
            int d = k - 1;
            for (; d >= 0; d--)
            {
              if (++choice[d] < units[order[combination[d]]].choices.size()) break;
              choice[d] = 0;
            }
            if (d < 0) break;
          }
          int i = k - 1;
          while (i >= 0 && combination[i] == units.size() - k + i) i--;
          if (i < 0) break;
          combination[i]++;
          for (int j = i + 1; j < k; j++) combination[j] = combination[j - 1] + 1;
        }
      }
      pp.form_entries.push_back(std::move(form_entries));
    }
    return pp;
  };

  // the decoy of a target: one permutation of its residues within the segments between form boundaries, keeping the first and the last
  // residue of every form and every non-standard residue; dst: target position -> decoy position (0-based)
  auto makeDecoy = [&](const ParsedProtein& target, std::vector<int>& dst) {
    ParsedProtein decoy = target;
    decoy.decoy = true;
    decoy.identifier = "DECOY_" + target.identifier;
    const int L = (int)target.sequence.size();
    std::vector<bool> fixed_point(L, false);
    std::set<int> cuts {0, L};
    for (const auto& f : target.forms)
    {
      fixed_point[f.start - 1] = true;
      fixed_point[f.end - 1] = true;
      cuts.insert(f.start - 1);
      cuts.insert(f.end);
    }
    for (int p = 0; p < L; p++)
    {
      if (movable_residues.find(target.sequence[p]) == std::string::npos) fixed_point[p] = true;
    }
    std::vector<int> src(L);
    for (int p = 0; p < L; p++) src[p] = p;
    std::mt19937_64 rng(fnv1a64(target.accession));
    for (auto it = cuts.begin(); it != cuts.end() && std::next(it) != cuts.end(); ++it)
    {
      const int a = *it, b = *std::next(it);
      std::vector<int> free_positions;
      for (int p = a; p < b; p++)
      {
        if (! fixed_point[p]) free_positions.push_back(p);
      }
      std::vector<int> values = free_positions;
      for (Size i = values.size(); i-- > 1;)
      {
        const Size j = (Size)(rng() % (std::uint64_t)(i + 1));
        std::swap(values[i], values[j]);
      }
      for (Size t = 0; t < free_positions.size(); t++) src[free_positions[t]] = values[t];
    }
    dst.assign(L, 0);
    for (int p = 0; p < L; p++)
    {
      decoy.sequence[p] = target.sequence[src[p]];
      dst[src[p]] = p;
    }
    for (auto& form_entries : decoy.form_entries)
    {
      for (auto& fixed : form_entries)
      {
        for (auto& [m, kind] : fixed)
        {
          m.position = dst[m.position - 1] + 1;
          if (m.bond > 0) m.bond = dst[m.bond - 1] + 1;
        }
        std::stable_sort(fixed.begin(), fixed.end(), [](const auto& a, const auto& b) { return a.first.position < b.first.position; });
      }
    }
    return decoy;
  };

  auto emit = [&](const ParsedProtein& pp) {
    const Size p = proteins.size();
    DatabaseProtein protein;
    protein.identifier = pp.identifier;
    protein.description = pp.description;
    protein.sequence = pp.sequence;
    protein.decoy = pp.decoy;
    for (const auto& f : pp.forms) protein.forms.emplace_back(f.start, f.end);
    for (Size fi = 0; fi < pp.forms.size(); fi++)
    {
      const auto& f = pp.forms[fi];
      const String letters = pp.sequence.substr(f.start - 1, f.end - f.start + 1);
      for (const auto& fixed : pp.form_entries[fi])
      {
        ProteoformEntry entry;
        entry.protein = p;
        entry.offset = f.start - 1;
        entry.length = f.end - f.start + 1;
        entry.id = pp.identifier + " " + String(f.start) + "-" + String(f.end);
        for (const auto& [m, kind] : fixed)
        {
          FixedModification fm;
          fm.position = m.position - f.start;
          fm.mass = m.mass;
          fm.accession = m.accession;
          fm.name = m.name;
          fm.kind = kind;
          fm.sites = m.sites;
          fm.bond = m.bond > 0 ? m.bond - f.start : -1;
          entry.fixed.push_back(fm);
        }
        try
        {
          ProForma pf(letters);
          for (const auto& fm : entry.fixed) pf.addModification(fm.position, fm.position, fm.mass, proformaTag(fm.accession, fm.mass));
          entry.proforma = pf.toString();
        }
        catch (const Exception::BaseException&)
        {
          entry.proforma = ""; // a residue ProForma has no letter for (e.g. X)
        }
        protein.entries.push_back(entries.size());
        fasta_entries.emplace_back(pp.identifier, pp.description, letters);
        entries.push_back(std::move(entry));
      }
    }
    proteins.push_back(std::move(protein));
  };

  // Every entry parsed first, so that the decoys of the file are paired by the identifiers as rebuilt (a PEFF sp:ACC with \ID is
  // sp|ACC|ID, as a FASTA decoy of it is named), as FLASHTnT pairs the entries of a FASTA file: a decoy with the first target of its name
  // without DECOY_ or DECOY, a target with the first decoy named DECOY_ or DECOY + its name. A decoy without such a target is paired with
  // the target of its accession (if one target has it), and that target with the first such decoy. Entries without residues are left out.
  std::vector<ParsedProtein> parsed;
  parsed.reserve(input.size());
  for (const auto& fe : input)
  {
    if (fe.sequence.empty())
    {
      statistics_.empty_sequences++;
      continue;
    }
    parsed.push_back(parseProtein(fe, fe.identifier.hasPrefix("DECOY")));
  }
  std::unordered_map<std::string, Size> first_of_name;
  std::unordered_map<std::string, long> target_of_accession; // -1: the accession of several targets
  for (Size i = 0; i < parsed.size(); i++)
  {
    first_of_name.emplace(parsed[i].identifier, i);
    if (parsed[i].decoy) continue;
    const auto [it, inserted] = target_of_accession.emplace(parsed[i].accession, (long)i);
    if (! inserted && it->second != (long)i) it->second = -1;
  }
  std::vector<long> partner_of(parsed.size(), -1); // indices into parsed
  for (Size i = 0; i < parsed.size(); i++)
  {
    if (! parsed[i].decoy) continue;
    statistics_.input_decoys++;
    if (const auto it = first_of_name.find(withoutDecoyPrefix(parsed[i].identifier)); it != first_of_name.end() && ! parsed[it->second].decoy)
    {
      partner_of[i] = (long)it->second;
    }
    else if (const auto ia = target_of_accession.find(parsed[i].accession); ia != target_of_accession.end() && ia->second >= 0)
    {
      partner_of[i] = ia->second;
    }
    else { statistics_.unpaired_decoys++; }
  }
  std::vector<long> decoy_by_accession(parsed.size(), -1); // per target: the first decoy paired with it by the accession
  for (Size i = 0; i < parsed.size(); i++)
  {
    if (! parsed[i].decoy || partner_of[i] < 0 || decoy_by_accession[(Size)partner_of[i]] >= 0) continue;
    if (withoutDecoyPrefix(parsed[i].identifier) != parsed[(Size)partner_of[i]].identifier) decoy_by_accession[(Size)partner_of[i]] = (long)i;
  }
  for (Size i = 0; i < parsed.size(); i++)
  {
    if (parsed[i].decoy) continue;
    for (const String& prefix : {"DECOY_", "DECOY"})
    {
      const auto it = first_of_name.find(prefix + parsed[i].identifier);
      if (it == first_of_name.end()) continue;
      partner_of[i] = (long)it->second;
      break;
    }
    if (partner_of[i] < 0) partner_of[i] = decoy_by_accession[i];
  }

  // the proteins in the file order; a target without a decoy in the file gets one, right after it
  std::vector<Size> protein_of_parsed(parsed.size());
  for (Size i = 0; i < parsed.size(); i++)
  {
    protein_of_parsed[i] = proteins.size();
    emit(parsed[i]);
    if (parsed[i].decoy) continue;
    statistics_.input_targets++;
    if (partner_of[i] >= 0) continue;
    std::vector<int> dst;
    emit(makeDecoy(parsed[i], dst));
    proteins[protein_of_parsed[i]].partner = (long)protein_of_parsed[i] + 1;
    proteins.back().partner = (long)protein_of_parsed[i];
    statistics_.generated_decoys++;
  }
  // the decoys of the file: partners, and whether their entries mirror the target's (the same number, the same masses), as the
  // target-decoy competition of the entries assumes
  auto entryMasses = [&](Size p) {
    std::vector<double> masses;
    for (Size e : proteins[p].entries)
    {
      double mass = 0;
      for (char c : fasta_entries[e].sequence) mass += residueMass(c);
      for (const auto& fm : entries[e].fixed) mass += fm.mass;
      masses.push_back(mass);
    }
    std::sort(masses.begin(), masses.end());
    return masses;
  };
  for (Size i = 0; i < parsed.size(); i++)
  {
    if (partner_of[i] < 0) continue;
    proteins[protein_of_parsed[i]].partner = (long)protein_of_parsed[(Size)partner_of[i]];
    if (! parsed[i].decoy) continue;
    const auto target_masses = entryMasses(protein_of_parsed[(Size)partner_of[i]]), decoy_masses = entryMasses(protein_of_parsed[i]);
    bool mirrored = target_masses.size() == decoy_masses.size();
    for (Size k = 0; mirrored && k < target_masses.size(); k++) mirrored = std::abs(target_masses[k] - decoy_masses[k]) <= 1e-3;
    statistics_.unmirrored_decoys += ! mirrored;
  }
  finish_();
}
} // namespace OpenMS
