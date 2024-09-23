// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong, Ayesha Feroz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/FORMAT/FLASHTnTFile.h>
#include <OpenMS/CHEMISTRY/ProForma.h>  // added for ProForma

namespace OpenMS
{
/**
  @brief FLASHTnT output *.tsv file format
   @ingroup FileIO
**/
void FLASHTnTFile::writeTagHeader(std::fstream& fs)
{
  fs << "TagIndex\tScan\tRetentionTime\tProteoformIndex\tProteinAccession\tProteinDescription\tTagSequence\tNmass\tCmass\tStartPosition\tDeltaMass\tL"
        "ength\tDeNovoScore\tMasses\tMassScores\n";
}

/// write header line for PrSM file
void FLASHTnTFile::writePrSMHeader(std::fstream& fs)
{
  fs << "ProteoformIndex\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tProteoformMass\tDatabaseSequence\tProteinSequence\tProf"
        "orma\tMatchingFragments\tCoverage(%)\tStartPosition\tEndPosition"
        "\tTagCount\tTagIndices\tModCount\tModMass\tModID\tModAccession\tModStart\tModEnd\tScore\tPrSMLevelQvalue\tProteoformLevelQvalue\n";
}

/// write header line for Proteoform file
void FLASHTnTFile::writeProHeader(std::fstream& fs)
{
  fs << "ProteoformIndex\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tProteoformMass\tDatabaseSequence\tProteinSequence\tProf"
        "orma\tMatchingFragments\tCoverage(%)\tStartPosition\tEndPosition"
        "\tTagCount\tTagIndices\tModCount\tModMass\tModID\tModAccession\tModStart\tModEnd\tScore\tProteoformLevelQvalue\n";
}

/// write the features in regular file output
void FLASHTnTFile::writeTags(const FLASHTnTAlgorithm& tnt, double flanking_mass_tol, std::fstream& fs)
{
  auto tags = std::vector<FLASHHelperClasses::Tag>();
  tnt.getTags(tags);

  for (int c = 0; c < 2; c++)
  {
    for (const auto& tag : tags)
    {
      auto hits = std::vector<ProteinHit>();
      tnt.getProteoformHitsMatchedBy(tag, hits);
      if (c == 0 && hits.empty()) continue;
      if (c == 1 && ! hits.empty()) continue;

      String acc = "";
      String description = "";
      String hitindices = "";
      String positions = "";
      String delta_masses = "";
      for (const auto& hit : hits)
      {
        if (! acc.empty()) acc += ";";
        if (! description.empty()) description += ";";
        if (! hitindices.empty()) hitindices += ";";
        if (! positions.empty()) positions += ";";
        if (! delta_masses.empty()) delta_masses += ";";

        acc += hit.getAccession();

        String proteindescription = hit.getDescription();
        if (proteindescription.empty()) { proteindescription = " "; }
        description += proteindescription;
        hitindices += (String)hit.getMetaValue("Index");

        auto pos = std::vector<int>();
        auto masses = std::vector<double>();
        FLASHTaggerAlgorithm::getMatchedPositionsAndFlankingMassDiffs(pos, masses, flanking_mass_tol, hit, tag);
        if (pos.size() != 0) { positions += std::to_string(pos[0]); }
        if (masses.size() != 0) { delta_masses += std::to_string(masses[0]); }
      }

      fs << tag.getIndex() << "\t" << tag.getScan() << "\t" << tag.getRetentionTime() << "\t" << hitindices << "\t" << acc << "\t" << description
         << "\t" << tag.getSequence() << "\t" << std::to_string(tag.getNtermMass()) << "\t" << std::to_string(tag.getCtermMass()) << "\t" << positions
         << "\t" << delta_masses << "\t" << tag.getLength() << "\t" << tag.getScore() << "\t";

      for (const auto& mz : tag.getMzs())
      {
        fs << std::to_string(mz) << ",";
      }
      fs << "\t";
      for (size_t i = 0; i <= tag.getLength(); i++) // Fixed signed/unsigned comparison issue
      {
        fs << std::to_string(tag.getScore(i)) << ",";
      }
      fs << "\n";
    }
  }
}

String generateProFormaString(const String& sequence,
                              const std::vector<double>& mod_masses,
                              const std::vector<int>& mod_starts,
                              const std::vector<int>& mod_ends,
                              const std::vector<String>& mod_ids,
                              const std::vector<String>& mod_accs)
{
  ProForma proforma(AASequence::fromString(sequence)); // Create ProForma object

  // Map to hold ranges that require parentheses
  std::map<std::pair<int, int>, String> range_modifications;

  // Loop through modifications and add them to the sequence
  for (size_t i = 0; i < mod_masses.size(); ++i)
  {
    // Check if the modification applies to a range
    if (mod_starts[i] != mod_ends[i])
    {
      // Mark range modifications using round brackets
      String range_mod = "(" + sequence.substr(mod_starts[i], mod_ends[i] - mod_starts[i]) + ")[";

      // Add modification information
      if (! mod_ids[i].empty()) { range_mod += mod_ids[i]; }
      else if (mod_masses[i] != 0.0) { range_mod += (mod_masses[i] > 0 ? "+" : "") + std::to_string(mod_masses[i]); }

      if (! mod_accs[i].empty()) { range_mod += ":" + mod_accs[i]; }

      range_mod += "]";
      range_modifications[{mod_starts[i] - 1, mod_ends[i] - 1}] = range_mod;
    }
    else
    {
      // For single position modifications, add them directly
      proforma.addModification(mod_starts[i] - 1, mod_ids[i], mod_masses[i]);
    }
  }

  // Build the ProForma string with range modifications
  String proforma_str = proforma.toProFormaString();

  // Insert range modifications into the ProForma string
  for (const auto& range_mod : range_modifications)
  {
    // Insert the range modification at the correct position
    int start_pos = range_mod.first.first; // Start of the range
    int end_pos = range_mod.first.second;  // End of the range

    // Find where to insert the range modification in the ProForma string
    proforma_str.insert(start_pos, range_mod.second);
  }

  return proforma_str; // Return the final ProForma string
}




void OpenMS::FLASHTnTFile::writePrSMs(const std::vector<ProteinHit>& hits, std::fstream& fs)
{
  for (const auto& hit : hits)
  {
    if (! hit.metaValueExists("Index")) continue;
    String tagindices = "", modmasses = "", modstarts = "", modends = "", modids = "", modaccs = "";

    int cntr = 0;
    std::vector<FLASHHelperClasses::Tag> tags;
    std::vector<int> indices;
    if (hit.metaValueExists("TagIndices")) indices = (std::vector<int>)hit.getMetaValue("TagIndices").toIntList();
    for (const int& index : indices)
    {
      if (! tagindices.empty()) tagindices += ";";
      tagindices += std::to_string(index);
      cntr++;
    }

    std::vector<double> mod_masses = hit.getMetaValue("Modifications");
    std::vector<int> mod_starts = hit.getMetaValue("ModificationStarts");
    std::vector<int> mod_ends = hit.getMetaValue("ModificationEnds");
    std::vector<String> mod_ids = hit.getMetaValue("ModificationIDs");
    std::vector<String> mod_accs = hit.getMetaValue("ModificationACCs");

    for (size_t i = 0; i < mod_masses.size(); i++) // Fixed signed/unsigned comparison issue
    {
      if (! modmasses.empty()) modmasses += ";";
      modmasses += std::to_string(mod_masses[i]);

      if (! modstarts.empty()) modstarts += ";";
      modstarts += std::to_string(mod_starts[i] + 1);

      if (! modends.empty()) modends += ";";
      modends += std::to_string(mod_ends[i] + 1);

      if (! modids.empty()) modids += ";";
      modids += mod_ids[i];

      if (! modaccs.empty()) modaccs += ";";
      modaccs += mod_accs[i];
    }

    int start = hit.getMetaValue("StartPosition");
    int end = hit.getMetaValue("EndPosition");

    int start_in_seq = start < 0 ? 0 : (start - 1);
    int end_in_seq = end < 0 ? hit.getSequence().size() : end;

    // Use ProForma for sequence generation
    String proformaStr = "";//generateProFormaString(hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq), mod_masses, mod_starts, mod_ends,
                                                //mod_ids, mod_accs);

    fs << hit.getMetaValue("Index") << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getMetaValue("RT") << "\t" << hit.getMetaValue("NumMass")
       << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getMetaValue("Mass") << "\t" << hit.getSequence() << "\t"
       << hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << hit.getMetaValue("MatchedAA") << "\t"
       << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << cntr << "\t" << tagindices << "\t" << mod_masses.size() << "\t"
       << modmasses << "\t" << modids << "\t" << modaccs << "\t" << modstarts << "\t" << modends << "\t" << hit.getScore() << "\t"
       << std::to_string((hit.metaValueExists("qvalue") ? (double)hit.getMetaValue("qvalue") : -1)) << "\t"
       << std::to_string((hit.metaValueExists("proqvalue") ? (double)hit.getMetaValue("proqvalue") : -1)) << "\n";
  }
}

void OpenMS::FLASHTnTFile::writeProteoforms(const std::vector<ProteinHit>& hits, std::fstream& fs, double pro_fdr)
{
  for (const auto& hit : hits)
  {
    if (! hit.metaValueExists("Index")) continue;
    if (! hit.metaValueExists("Representative")) continue;
    if (hit.metaValueExists("proqvalue") && (double)hit.getMetaValue("proqvalue") > pro_fdr) continue;
    String tagindices = "", modmasses = "", modstarts = "", modends = "", modids = "", modaccs = "";

    int cntr = 0;
    std::vector<FLASHHelperClasses::Tag> tags;
    std::vector<int> indices;
    if (hit.metaValueExists("TagIndices")) indices = (std::vector<int>)hit.getMetaValue("TagIndices").toIntList();
    for (const int& index : indices)
    {
      if (! tagindices.empty()) tagindices += ";";
      tagindices += std::to_string(index);
      cntr++;
    }

    std::vector<double> mod_masses = hit.getMetaValue("Modifications");
    std::vector<int> mod_starts = hit.getMetaValue("ModificationStarts");
    std::vector<int> mod_ends = hit.getMetaValue("ModificationEnds");
    std::vector<String> mod_ids = hit.getMetaValue("ModificationIDs");
    std::vector<String> mod_accs = hit.getMetaValue("ModificationACCs");

    for (size_t i = 0; i < mod_masses.size(); i++) // Fixed signed/unsigned comparison issue
    {
      if (! modmasses.empty()) modmasses += ";";
      modmasses += std::to_string(mod_masses[i]);

      if (! modstarts.empty()) modstarts += ";";
      modstarts += std::to_string(mod_starts[i] + 1);

      if (! modends.empty()) modends += ";";
      modends += std::to_string(mod_ends[i] + 1);

      if (! modids.empty()) modids += ";";
      modids += mod_ids[i];

      if (! modaccs.empty()) modaccs += ";";
      modaccs += mod_accs[i];
    }

    int start = hit.getMetaValue("StartPosition");
    int end = hit.getMetaValue("EndPosition");

    int start_in_seq = start < 0 ? 0 : (start - 1);
    int end_in_seq = end < 0 ? hit.getSequence().size() : end;

    // Use ProForma to generate ProForma string
    String proformaStr = "";//generateProFormaString(hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq), mod_masses, mod_starts, mod_ends,
                                                //mod_ids, mod_accs);

    fs << hit.getMetaValue("Index") << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getMetaValue("RT") << "\t" << hit.getMetaValue("NumMass")
       << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getMetaValue("Mass") << "\t" << hit.getSequence() << "\t"
       << hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << hit.getMetaValue("MatchedAA") << "\t"
       << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << cntr << "\t" << tagindices << "\t" << mod_masses.size() << "\t"
       << modmasses << "\t" << modids << "\t" << modaccs << "\t" << modstarts << "\t" << modends << "\t" << hit.getScore() << "\t"
       << std::to_string((hit.metaValueExists("proqvalue") ? (double)hit.getMetaValue("proqvalue") : -1)) << "\n";
  }
}

}
// namespace OpenMS