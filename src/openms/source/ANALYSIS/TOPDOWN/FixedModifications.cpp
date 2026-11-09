// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FixedModifications.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/ProteoformDatabase.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <cctype>

namespace OpenMS
{
  const FixedModifications& FixedModifications::none()
  {
    static const FixedModifications empty_table;
    return empty_table;
  }

  FixedModifications FixedModifications::fromNames(const std::vector<std::string>& names)
  {
    FixedModifications table;
    for (const std::string& given : names)
    {
      const String name = String(given).trim();
      if (name.empty()) continue;
      const ResidueModification* mod = nullptr;
      try
      {
        mod = FLASHHelperClasses::resolveModification(name); // a name or a Unimod accession, with its residue
      }
      catch (const Exception::BaseException& e)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "fixed_mod " + name + ": " + e.what(), name);
      }
      const char origin = mod->getOrigin();
      if (mod->getTermSpecificity() != ResidueModification::ANYWHERE || ! std::isupper((unsigned char)origin) || origin == 'X')
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "fixed_mod " + name + ": only a modification of a residue anywhere can be fixed (not a terminal one).", name);
      }
      const unsigned char u = (unsigned char)origin;
      if (table.fixed_[u])
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "fixed_mod: two fixed modifications of residue " + String(origin) + ".", name);
      }
      // the mass the FLASHTnT tool's ResidueDB overwrite gave the residue: its full weight plus the modification
      Residue residue(*ResidueDB::getInstance()->getResidue(u));
      residue.setMonoWeight(residue.getMonoWeight(Residue::Full) + mod->getDiffMonoMass());
      table.fixed_[u] = true;
      table.delta_[u] = mod->getDiffMonoMass();
      table.modified_[u] = residue.getMonoWeight(Residue::Internal);
      Entry entry;
      entry.residue = origin;
      entry.mass = mod->getDiffMonoMass();
      entry.accession = mod->getUniModAccession();
      entry.accession.substitute("UniMod:", "UNIMOD:");
      entry.name = name;
      table.entries_.push_back(entry);
    }
    return table;
  }

  FixedModifications FixedModifications::fromProFormaPrefix(const String& proforma)
  {
    // <[tag]@residues> groups before the sequence; tag = "+mass|UNIMOD:n" (ProteoformDatabase::proformaTag)
    std::vector<std::string> names;
    Size pos = 0;
    while (pos < proforma.size() && proforma[pos] == '<')
    {
      const Size close = proforma.find('>', pos);
      const Size tag_end = proforma.find("]@", pos);
      if (close == String::npos || tag_end == String::npos || tag_end > close || pos + 1 >= proforma.size() || proforma[pos + 1] != '[')
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "a malformed ProForma global modification", proforma);
      }
      const String tag = proforma.substr(pos + 2, tag_end - pos - 2);
      const String residues = proforma.substr(tag_end + 2, close - tag_end - 2);
      String accession;
      std::vector<String> parts;
      tag.split('|', parts);
      if (parts.empty()) parts.push_back(tag);
      for (const String& part : parts)
      {
        if (String(part).toUpper().hasPrefix("UNIMOD:")) accession = part;
      }
      if (accession.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "a fixed modification without a Unimod accession: <[" + tag + "]@" + residues + ">", proforma);
      }
      std::vector<String> sites;
      residues.split(',', sites);
      if (sites.empty()) sites.push_back(residues);
      for (const String& site : sites) names.push_back(accession + " (" + String(site).trim() + ")");
      pos = close + 1;
    }
    return fromNames(names);
  }

  bool FixedModifications::empty() const
  {
    return entries_.empty();
  }

  double FixedModifications::mass(char c) const
  {
    const unsigned char u = (unsigned char)c;
    if (fixed_[u]) return modified_[u];
    return ResidueDB::getInstance()->getResidue(u)->getMonoWeight(Residue::Internal);
  }

  double FixedModifications::delta(char c) const
  {
    return delta_[(unsigned char)c];
  }

  double FixedModifications::internalMass(const AASequence& seq) const
  {
    double m = seq.getMonoWeight(Residue::Internal);
    if (entries_.empty()) return m;
    for (const Residue& r : seq)
    {
      const String& code = r.getOneLetterCode();
      if (! code.empty()) m += delta(code[0]);
    }
    return m;
  }

  const std::vector<FixedModifications::Entry>& FixedModifications::entries() const
  {
    return entries_;
  }

  std::vector<std::string> FixedModifications::names() const
  {
    std::vector<std::string> out;
    out.reserve(entries_.size());
    for (const Entry& e : entries_) out.push_back(e.name);
    return out;
  }

  String FixedModifications::proformaPrefix() const
  {
    String prefix;
    for (const Entry& e : entries_) prefix += "<[" + ProteoformDatabase::proformaTag(e.accession, e.mass) + "]@" + String(e.residue) + ">";
    return prefix;
  }

  bool FixedModifications::operator==(const FixedModifications& other) const
  {
    return fixed_ == other.fixed_ && delta_ == other.delta_;
  }

  bool FixedModifications::operator!=(const FixedModifications& other) const
  {
    return ! (*this == other);
  }
} // namespace OpenMS
