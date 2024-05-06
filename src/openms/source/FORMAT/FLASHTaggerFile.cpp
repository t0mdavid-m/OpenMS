// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/FORMAT/FLASHTaggerFile.h>

namespace OpenMS
{
  /**
    @brief FLASHDeconv Tagger output *.tsv file format
     @ingroup FileIO
  **/
  void FLASHTaggerFile::writeTagHeader(std::fstream& fs)
  {
    fs << "Scan\ttagIndex\tProteinIndex\tProteinAccession\tProteinDescription\tTagSequence\tNmass\tCmass\tStartPos\tDeltaMass\tLength\tDeNovoScore\tMasses\tMassScores\t"
          "\n";
  }

  /// write header line for topFD feature file
  void FLASHTaggerFile::writeProteinHeader(std::fstream& fs)
  {
    fs << "Scan\tProteinIndex\tProteinAccession\tProteinDescription\tProteinSequence\tMatchedAminoAcidCount\tCoverage(%)"
          "\tProteinScore\tProteinQvalue\tTagCount\tTagIndices\n";
  }

  /// write the features in regular file output
  void FLASHTaggerFile::writeTags(const FLASHTaggerAlgorithm& tagger, std::fstream& fs)
  {
    for (int c = 0; c < 2; c++)
    {
      for (int scan : tagger.getScans())
      {
        for (const auto& tag : tagger.getTagsAt(scan))
        {
          auto hits = tagger.getProteinHitsMatchedBy(tag);

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
            description += hit.getDescription();
            hitindices += std::to_string(tagger.getProteinIndex(hit, tag.getScan()));

            auto seqposition = tagger.getMatchedPositions(hit, tag);
            if (seqposition.size() != 0) { positions += std::to_string(seqposition[0]); }

            auto delta_mass = tagger.getDeltaMasses(hit, tag);
            if (delta_mass.size() != 0) { delta_masses += std::to_string(delta_mass[0]); }
          }

          fs << tag.getScan() << "\t" << tagger.getTagIndex(tag) << "\t" << hitindices << "\t" << acc << "\t" << description << "\t"
             << tag.getSequence() << "\t" << std::to_string(tag.getNtermMass()) << "\t" << std::to_string(tag.getCtermMass()) << "\t" << positions
             << "\t" << delta_masses << "\t" << tag.getLength() << "\t" << tag.getScore() << "\t";

          for (const auto& mz : tag.getMzs())
          {
            fs << std::to_string(mz) << ",";
          }
          fs << "\t";
          for (auto i = 0; i < tag.getLength(); i++)
          {
            fs << std::to_string(tag.getScore(i)) << ",";
          }
          fs << "\n";
        }
      }
    }
  }

  void FLASHTaggerFile::writeProteins(const FLASHTaggerAlgorithm& tagger, std::fstream& fs)
  {
    for (const int scan : tagger.getScans())
    {
      for (const auto& hit : tagger.getProteinHitsAt(scan))
      {
        String tagindices = "";
        int cntr = 0;
        for (const auto& tag : tagger.getTagsMatchingTo(hit, scan))
        {
          if (! tagindices.empty()) tagindices += ";";
          tagindices += std::to_string(tagger.getTagIndex(tag));
          cntr++;
        }

        fs << scan << "\t" << tagger.getProteinIndex(hit, scan) << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getSequence() << "\t"
           << hit.getMetaValue("MatchedAA") << "\t" << 100.0 * hit.getCoverage() << "\t" << hit.getScore() << "\t"
           << std::to_string((hit.metaValueExists("qvalue")? (double)hit.getMetaValue("qvalue") : 0)) << "\t" << cntr << "\t" << tagindices << "\n";
      }
    }
  }
} // namespace OpenMS