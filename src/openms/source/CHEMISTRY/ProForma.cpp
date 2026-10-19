// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ProForma.h>

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>

namespace OpenMS
{
namespace
{
  const std::string residue_letters = "ACDEFGHIKLMNPQRSTVWYUO";

  [[noreturn]] void fail(const String& input, const String& message)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input, message);
  }

  /// the content of the bracket opening at s[pos]; pos moves past its closing bracket. Brackets of the same kind nest (Cation:Ca[II]).
  String readBracket(const String& s, Size& pos, char open, char close)
  {
    int depth = 0;
    for (Size i = pos; i < s.size(); ++i)
    {
      if (s[i] == open) depth++;
      else if (s[i] == close && --depth == 0)
      {
        String content = s.substr(pos + 1, i - pos - 1);
        pos = i + 1;
        return content;
      }
    }
    fail(s, String("unclosed '") + open + "'");
  }

  /// the monoisotopic masses of the ProForma glycan names
  double glycanMass(const String& composition)
  {
    // longest names first, as several share a prefix (HexNAc, HexN, HexA, HexS, Hex)
    static const std::vector<std::pair<String, String>> monosaccharides = {
      {"HexNAcS", "C8H13NO8S"}, {"HexNAc", "C8H13NO5"}, {"NeuAc", "C11H17NO8"}, {"NeuGc", "C11H17NO9"}, {"HexS", "C6H10O8S"},
      {"HexP", "C6H11O8P"},     {"HexN", "C6H11NO4"},   {"HexA", "C6H8O6"},     {"dHex", "C6H10O4"},   {"Hex", "C6H10O5"},
      {"Fuc", "C6H10O4"},       {"Neu", "C9H15NO7"},    {"Pen", "C5H8O4"},      {"Kdn", "C9H14O8"},    {"Phosphate", "HPO3"},
      {"Sulfate", "SO3"}};
    double mass = 0;
    Size pos = 0;
    while (pos < composition.size())
    {
      if (composition[pos] == ' ')
      {
        pos++;
        continue;
      }
      bool found = false;
      for (const auto& [name, formula] : monosaccharides)
      {
        if (composition.compare(pos, name.size(), name) != 0) continue;
        pos += name.size();
        Size digits = pos;
        while (digits < composition.size() && std::isdigit((unsigned char)composition[digits])) digits++;
        const int count = digits > pos ? std::stoi(composition.substr(pos, digits - pos)) : 1;
        pos = digits;
        mass += count * EmpiricalFormula(formula).getMonoWeight();
        found = true;
        break;
      }
      if (! found) fail(composition, "unknown monosaccharide in glycan composition");
    }
    return mass;
  }

  /// a modification name or accession looked up in ModificationsDB (residue-specific first)
  bool lookUpMass(const String& name, char residue, double& mass)
  {
    const auto* db = ModificationsDB::getInstance();
    if (! db->has(name)) return false; // searchModifications would log every miss
    std::set<const ResidueModification*> mods;
    if (residue != 0) db->searchModifications(mods, name, String(residue));
    if (mods.empty()) db->searchModifications(mods, name);
    if (mods.empty()) return false;
    mass = (*mods.begin())->getDiffMonoMass();
    return true;
  }

  /// a number with an optional sign; false if the whole string is not a number
  bool toNumber(const String& s, double& value)
  {
    if (s.empty()) return false;
    try
    {
      Size used = 0;
      value = std::stod(s, &used);
      return used == s.size();
    }
    catch (const std::exception&)
    {
      return false;
    }
  }

  /// the mass of one alternative of a tag; false if it gives none (INFO, an empty localisation group reference)
  bool alternativeMass(String part, char residue, double& mass)
  {
    part.trim();
    // a localisation group or a score: Phospho#g1(0.9) -> Phospho
    const Size hash = part.find('#');
    if (hash != String::npos)
    {
      const String suffix = part.substr(hash + 1);
      if (suffix.hasPrefix("XL") || suffix.hasPrefix("BRANCH")) fail(part, "cross-links and branches are not supported");
      part = part.substr(0, hash);
    }
    else if (! part.empty() && part.back() == ')')
    {
      const Size open = part.rfind('(');
      double score;
      if (open != String::npos && toNumber(part.substr(open + 1, part.size() - open - 2), score)) part = part.substr(0, open);
    }
    part.trim();
    if (part.empty()) return false;
    const Size colon = part.find(':');
    const String prefix = colon == String::npos ? String() : String(part.substr(0, colon)).toUpper();
    const String value = colon == String::npos ? part : String(part.substr(colon + 1));
    if (prefix == "INFO") return false;
    if (toNumber(part, mass)) return true;
    if (prefix == "OBS")
    {
      if (! toNumber(value, mass)) fail(part, "observed mass is not a number");
      return true;
    }
    if (prefix == "FORMULA")
    {
      String formula = value;
      formula.substitute(" ", "");
      try
      {
        mass = EmpiricalFormula(formula).getMonoWeight();
      }
      catch (const Exception::BaseException&)
      {
        fail(part, "invalid formula");
      }
      return true;
    }
    if (prefix == "GLYCAN")
    {
      mass = glycanMass(value);
      return true;
    }
    if (prefix == "UNIMOD") return lookUpMass("UniMod:" + value, residue, mass);
    if (prefix == "MOD") return lookUpMass("MOD:" + value, residue, mass);
    if (prefix == "U" || prefix == "M") return lookUpMass(value, residue, mass);
    if (prefix == "R" || prefix == "RESID" || prefix == "X" || prefix == "XLMOD" || prefix == "G" || prefix == "GNO")
    {
      return false; // RESID, XL-MOD and GNOme are not in ModificationsDB: another alternative must give the mass
    }
    // a name without source, or one with a colon of its own (Cation:Na): Unimod, then PSI-MOD
    return lookUpMass(part, residue, mass);
  }

  /// the residues of a global fixed modification: <[tag]@C,M> or <[tag]@N-term:C> (only plain residues are supported)
  std::vector<char> globalResidues(const String& spec, const String& input)
  {
    std::vector<String> parts;
    spec.split(',', parts);
    std::vector<char> residues;
    for (auto p : parts)
    {
      p.trim();
      if (p.size() != 1 || residue_letters.find(p[0]) == std::string::npos) fail(input, "global modification on '" + p + "' is not supported");
      residues.push_back(p[0]);
    }
    if (residues.empty()) fail(input, "global modification without residues");
    return residues;
  }
} // namespace

ProForma::ProForma(const String& sequence) : sequence_(sequence)
{
  for (char c : sequence_)
  {
    if (residue_letters.find(c) == std::string::npos) fail(sequence, String("unsupported residue '") + c + "'");
  }
}

double ProForma::getTagMass(const String& tag, char residue)
{
  std::vector<String> alternatives;
  tag.split('|', alternatives);
  if (alternatives.empty()) alternatives.push_back(tag);
  for (const auto& alternative : alternatives)
  {
    double mass = 0;
    if (alternativeMass(alternative, residue, mass)) return mass;
  }
  fail(tag, "no mass for the modification");
}

String ProForma::formatMass(double mass)
{
  char buffer[64];
  std::snprintf(buffer, sizeof(buffer), "%+.6f", mass);
  String s(buffer);
  while (s.size() > 3 && s.back() == '0' && s[s.size() - 2] != '.') s.pop_back();
  return s;
}

void ProForma::addModification(int start, int end, double mass, const String& tag)
{
  if (start > end) std::swap(start, end);
  if (start < 0 || end >= (int)sequence_.size())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "modification outside the sequence", String(start) + "-" + String(end));
  }
  if (end > start)
  {
    for (const auto& m : mods_)
    {
      if (m.site != Site::RESIDUE || m.end == m.start || (m.start == start && m.end == end)) continue;
      // ranges are either the same or disjoint, and a range may hold no other range
      if (m.start <= end && start <= m.end)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "overlapping ranges", String(start) + "-" + String(end));
      }
    }
  }
  Modification m;
  m.site = Site::RESIDUE;
  m.start = start;
  m.end = end;
  m.mass = mass;
  m.tag = tag;
  mods_.push_back(m);
}

void ProForma::addModification(Site site, double mass, const String& tag)
{
  if (site == Site::RESIDUE) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "a residue modification needs its residues", "");
  Modification m;
  m.site = site;
  m.mass = mass;
  m.tag = tag;
  mods_.push_back(m);
}

double ProForma::getMonoMass() const
{
  double mass = EmpiricalFormula("H2O").getMonoWeight();
  for (char c : sequence_) mass += ResidueDB::getInstance()->getResidue(c)->getMonoWeight(Residue::Internal);
  for (const auto& m : mods_) mass += m.mass;
  return mass;
}

ProForma ProForma::fromString(const String& input)
{
  String s = input;
  s.trim();
  ProForma pf;
  Size pos = 0;

  // before the sequence: global, unknown-position, labile and N-terminal modifications, in any order
  std::vector<std::pair<String, std::vector<char>>> globals;
  std::vector<String> unknown, labile, n_term;
  while (pos < s.size())
  {
    if (s[pos] == '<')
    {
      const String content = readBracket(s, pos, '<', '>');
      if (content.empty() || content[0] != '[') fail(input, "isotope labels <...> are not supported");
      Size p = 0;
      const String tag = readBracket(content, p, '[', ']');
      if (p >= content.size() || content[p] != '@') fail(input, "global modification without '@'");
      globals.emplace_back(tag, globalResidues(content.substr(p + 1), input));
    }
    else if (s[pos] == '{') { labile.push_back(readBracket(s, pos, '{', '}')); }
    else if (s[pos] == '[')
    {
      std::vector<String> tags;
      while (pos < s.size() && s[pos] == '[') tags.push_back(readBracket(s, pos, '[', ']'));
      int count = 1;
      if (pos < s.size() && s[pos] == '^')
      {
        Size digits = ++pos;
        while (digits < s.size() && std::isdigit((unsigned char)s[digits])) digits++;
        if (digits == pos) fail(input, "'^' without a count");
        count = std::stoi(s.substr(pos, digits - pos));
        pos = digits;
        if (pos >= s.size() || s[pos] != '?') fail(input, "'^n' must be followed by '?'");
      }
      if (pos < s.size() && s[pos] == '?')
      {
        pos++;
        for (int c = 0; c < count; c++) unknown.insert(unknown.end(), tags.begin(), tags.end());
      }
      else if (pos < s.size() && s[pos] == '-')
      {
        pos++;
        n_term.insert(n_term.end(), tags.begin(), tags.end());
      }
      else fail(input, "a modification before the sequence must be followed by '?' or '-'");
    }
    else break;
  }

  // the sequence with its residue and range modifications, then the C-terminal modifications and the charge
  struct Pending
  {
    int start, end;
    String tag;
  };
  std::vector<Pending> residue_mods;
  std::vector<String> c_term;
  int range_start = -1;
  bool c_term_done = false;
  while (pos < s.size())
  {
    const char c = s[pos];
    if (c_term_done && c != '/') fail(input, "unexpected text after the C-terminal modifications");
    if (std::isupper((unsigned char)c))
    {
      if (residue_letters.find(c) == std::string::npos) fail(input, String("unsupported residue '") + c + "'");
      pf.sequence_ += c;
      pos++;
      const int index = (int)pf.sequence_.size() - 1;
      while (pos < s.size() && s[pos] == '[') residue_mods.push_back({index, index, readBracket(s, pos, '[', ']')});
    }
    else if (c == '(')
    {
      if (range_start >= 0) fail(input, "nested ranges are not supported");
      range_start = (int)pf.sequence_.size();
      pos++;
    }
    else if (c == ')')
    {
      if (range_start < 0 || range_start >= (int)pf.sequence_.size()) fail(input, "unmatched or empty range");
      pos++;
      if (pos >= s.size() || s[pos] != '[') fail(input, "a range must be followed by its modification");
      while (pos < s.size() && s[pos] == '[') residue_mods.push_back({range_start, (int)pf.sequence_.size() - 1, readBracket(s, pos, '[', ']')});
      range_start = -1;
    }
    else if (c == '-')
    {
      pos++;
      if (pos >= s.size() || s[pos] != '[') fail(input, "'-' must be followed by a C-terminal modification");
      while (pos < s.size() && s[pos] == '[') c_term.push_back(readBracket(s, pos, '[', ']'));
      c_term_done = true;
    }
    else if (c == '/')
    {
      pos++;
      Size digits = pos;
      if (digits < s.size() && (s[digits] == '+' || s[digits] == '-')) digits++;
      while (digits < s.size() && std::isdigit((unsigned char)s[digits])) digits++;
      if (digits == pos || digits != s.size()) fail(input, "only an integer charge is supported after '/'");
      pf.charge_ = std::stoi(s.substr(pos, digits - pos));
      pos = digits;
    }
    else if (c == '+') { fail(input, "chimeric spectra (+) are not supported"); }
    else fail(input, String("unexpected character '") + c + "'");
  }
  if (range_start >= 0) fail(input, "unclosed range");
  if (pf.sequence_.empty()) fail(input, "no sequence");

  // modifications; an empty localisation group reference ([#g1]) gives no mass and marks only another possible site
  auto has_mass = [&input](const String& tag) {
    if (tag.has('#'))
    {
      const String group = tag.substr(tag.find('#') + 1);
      if (group.hasPrefix("XL") || group.hasPrefix("BRANCH")) fail(input, "cross-links and branches are not supported");
    }
    std::vector<String> alternatives;
    tag.split('|', alternatives);
    for (auto a : alternatives)
    {
      a.trim();
      if (! a.empty() && a[0] != '#' && ! String(a).toUpper().hasPrefix("INFO:")) return true;
    }
    return false;
  };
  for (const auto& [tag, residues] : globals)
  {
    for (int i = 0; i < (int)pf.sequence_.size(); i++)
    {
      if (std::find(residues.begin(), residues.end(), pf.sequence_[i]) == residues.end()) continue;
      Modification m;
      m.start = m.end = i;
      m.mass = getTagMass(tag, pf.sequence_[i]);
      m.tag = tag;
      m.global = true;
      pf.mods_.push_back(m);
    }
  }
  for (const auto& tag : unknown) pf.addModification(Site::UNKNOWN, getTagMass(tag), tag);
  for (const auto& tag : labile) pf.addModification(Site::LABILE, getTagMass(tag), tag);
  for (const auto& tag : n_term) pf.addModification(Site::N_TERM, getTagMass(tag, pf.sequence_[0]), tag);
  for (const auto& r : residue_mods)
  {
    if (! has_mass(r.tag)) continue;
    pf.addModification(r.start, r.end, getTagMass(r.tag, r.start == r.end ? pf.sequence_[r.start] : 0), r.tag);
  }
  for (const auto& tag : c_term) pf.addModification(Site::C_TERM, getTagMass(tag, pf.sequence_.back()), tag);
  return pf;
}

String ProForma::toString() const
{
  auto tag_of = [](const Modification& m) { return "[" + (m.tag.empty() ? formatMass(m.mass) : m.tag) + "]"; };
  String out;

  // global fixed modifications, one per tag, with their residues
  std::map<String, std::set<char>> globals;
  for (const auto& m : mods_)
  {
    if (m.global) globals[tag_of(m)].insert(sequence_[m.start]);
  }
  for (const auto& [tag, residues] : globals)
  {
    out += "<" + tag + "@";
    bool first = true;
    for (char r : residues)
    {
      if (! first) out += ",";
      out += r;
      first = false;
    }
    out += ">";
  }
  bool any = false;
  for (const auto& m : mods_)
  {
    if (m.site == Site::UNKNOWN)
    {
      out += tag_of(m);
      any = true;
    }
  }
  if (any) out += "?";
  for (const auto& m : mods_)
  {
    if (m.site == Site::LABILE) out += "{" + (m.tag.empty() ? formatMass(m.mass) : m.tag) + "}";
  }
  any = false;
  for (const auto& m : mods_)
  {
    if (m.site == Site::N_TERM)
    {
      out += tag_of(m);
      any = true;
    }
  }
  if (any) out += "-";

  for (int i = 0; i < (int)sequence_.size(); i++)
  {
    int range_end = -1;
    for (const auto& m : mods_)
    {
      if (m.site == Site::RESIDUE && ! m.global && m.start == i && m.end > i) range_end = m.end;
    }
    if (range_end >= 0) out += "(";
    out += sequence_[i];
    for (const auto& m : mods_)
    {
      if (m.site == Site::RESIDUE && ! m.global && m.start == i && m.end == i) out += tag_of(m);
    }
    // close a range at its last residue, with all its modifications
    bool closed = false;
    for (const auto& m : mods_)
    {
      if (m.site != Site::RESIDUE || m.global || m.end != i || m.end == m.start) continue;
      if (! closed) out += ")";
      closed = true;
      out += tag_of(m);
    }
  }

  any = false;
  for (const auto& m : mods_)
  {
    if (m.site == Site::C_TERM)
    {
      if (! any) out += "-";
      out += tag_of(m);
      any = true;
    }
  }
  if (charge_ != 0) out += "/" + String(charge_);
  return out;
}
} // namespace OpenMS
