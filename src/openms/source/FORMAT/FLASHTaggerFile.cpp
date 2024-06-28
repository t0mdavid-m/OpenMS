// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/FORMAT/FLASHTaggerFile.h>

namespace OpenMS
{
  /**
    @brief FLASHDeconv Tagger output *.tsv file format
     @ingroup FileIO
  **/
  void FLASHTaggerFile::writeTagHeader(std::fstream& fs)
  {
    fs << "Scan\tTagIndex\tProteinIndex\tProteinAccession\tProteinDescription\tTagSequence\tNmass\tCmass\tStartPos\tDeltaMass\tLength\tDeNovoScore\tMas"
          "ses\tMassScores\t"
          "\n";
  }

  /// write header line for topFD feature file
  void FLASHTaggerFile::writeProteinHeader(std::fstream& fs)
  {
    fs << "ProteoformIndex\tScan\tProteinAccession\tProteinDescription\tProteinSequence\tMatchedAminoAcidCount\tCoverage(%)\tStartPosition\tEndPosition"
          "\tTagCount\tTagIndices\tModificationCount\tModifications\tModificationStarts\tModificationEnds\tProteoformScore\tProteinLevelQvalue\n";
  }

  /// write the features in regular file output
  void FLASHTaggerFile::writeTags(const FLASHTnTAlgorithm& tnt, double flanking_mass_tol, std::fstream& fs)
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
        if (c == 1 && !hits.empty()) continue;

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
          if (proteindescription.empty()) {
            proteindescription = " ";
          }
          description += proteindescription;
          hitindices += (String)hit.getMetaValue("Index");

          auto pos = std::vector<int>();
          auto masses = std::vector<double>();
          FLASHTaggerAlgorithm::getMatchedPositionsAndFlankingMassDiffs(pos, masses, flanking_mass_tol, hit, tag);
          if (pos.size() != 0) { positions += std::to_string(pos[0]); }
          if (masses.size() != 0) { delta_masses += std::to_string(masses[0]); }
        }

        fs << tag.getScan() << "\t" << tag.getIndex() << "\t" << hitindices << "\t" << acc << "\t" << description << "\t" << tag.getSequence() << "\t"
           << std::to_string(tag.getNtermMass()) << "\t" << std::to_string(tag.getCtermMass()) << "\t" << positions << "\t" << delta_masses << "\t"
           << tag.getLength() << "\t" << tag.getScore() << "\t";

        for (const auto& mz : tag.getMzs())
        {
          fs << std::to_string(mz) << ";";
        }
        fs << "\t";
        for (auto i = 0; i < tag.getLength(); i++)
        {
          fs << std::to_string(tag.getScore(i)) << ";";
        }
        fs << "\n";
      }
    }
  }

  void FLASHTaggerFile::writeProteins(const std::vector<ProteinHit>& hits, std::fstream& fs)
  {
    for (const auto& hit : hits)
    {
      if (!hit.metaValueExists("Index")) continue;
      String tagindices = "", modmasses = "", modstarts="", modeends="";

      int cntr = 0;
      std::vector<FLASHHelperClasses::Tag> tags;
      for (const auto& index : (std::vector<int>)hit.getMetaValue("TagIndices"))
      {
        if (! tagindices.empty()) tagindices += ";";
        tagindices += std::to_string(index);
        cntr++;
      }

      std::vector<double> mod_masses = hit.getMetaValue("Modifications");
      std::vector<int> mod_starts = hit.getMetaValue("ModificationStarts");
      std::vector<int> mod_ends = hit.getMetaValue("ModificationEnds");

      for (int i = 0; i < mod_masses.size(); i++)
      {
        if (! modmasses.empty()) modmasses += ";";
        modmasses += std::to_string(mod_masses[i]);

        if (! mod_starts.empty()) modstarts += ";";
        modstarts += std::to_string(mod_starts[i]);

        if (! mod_ends.empty()) modeends += ";";
        modeends += std::to_string(mod_ends[i]);
      }

      fs << hit.getMetaValue("Index") << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getSequence() << "\t"
         << hit.getMetaValue("MatchedAA") << "\t" << 100.0 * hit.getCoverage() << "\t" << hit.getMetaValue("StartPosition") <<"\t"<<hit.getMetaValue("EndPosition")  << "\t"
         <<  cntr << "\t" << tagindices << "\t"
         << mod_masses.size() <<"\t"<< modmasses <<"\t"<<modstarts <<"\t"  << modeends << "\t" << hit.getScore() << "\t" << std::to_string((hit.metaValueExists("qvalue") ? (double)hit.getMetaValue("qvalue") : -1)) << "\n";
      }
    }
  } // namespace OpenMS