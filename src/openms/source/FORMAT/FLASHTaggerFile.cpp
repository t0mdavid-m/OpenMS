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
    fs << "Scan\tTagIndex\tProteinIndex\tProteinAccession\tProteinDescription\tTagSequence\tNmass\tCmass\tStartPos\tDeltaMass\tLength\tDeNovoScore\tMas"
          "ses\tMassScores\t"
          "\n";
  }

  /// write header line for topFD feature file
  void FLASHTaggerFile::writeProteinHeader(std::fstream& fs)
  {
    fs << "ProteinIndex\tProteinAccession\tProteinDescription\tProteinSequence\tMatchedAminoAcidCount\tCoverage(%)"
          "\tProteinScore\tProteinQvalue\tTagCount\tTagIndices\tScans\n";
  }

  /// write the features in regular file output
  void FLASHTaggerFile::writeTags(const FLASHTaggerAlgorithm& tagger, std::fstream& fs)
  {
    for (int c = 0; c < 2; c++)
    {
      auto tags = std::vector<FLASHDeconvHelperStructs::Tag>();
      tagger.getTags(c == 0, tags);
      for (const auto& tag : tags)
      {
        auto hits = std::vector<ProteinHit>();
        tagger.getProteinHitsMatchedBy(tag, hits);

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
            proteindescription = " "
          }
          description += proteindescription;
          hitindices += std::to_string(tagger.getProteinIndex(hit));

          auto seqposition = tagger.getMatchedPositions(hit, tag);
          if (seqposition.size() != 0) { positions += std::to_string(seqposition[0]); }

          auto delta_mass = tagger.getDeltaMasses(hit, tag);
          if (delta_mass.size() != 0) { delta_masses += std::to_string(delta_mass[0]); }
        }

        fs << tag.getScan() << "\t" << tag.getIndex() << "\t" << hitindices << "\t" << acc << "\t" << description << "\t" << tag.getSequence() << "\t"
           << std::to_string(tag.getNtermMass()) << "\t" << std::to_string(tag.getCtermMass()) << "\t" << positions << "\t" << delta_masses << "\t"
           << tag.getLength() << "\t" << tag.getScore() << "\t";

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

  void FLASHTaggerFile::writeProteins(const FLASHTaggerAlgorithm& tagger, std::fstream& fs)
  {
    auto hits = std::vector<ProteinHit>();
    tagger.getProteinHits(hits);
    for (const auto& hit : hits)
    {
      String tagindices = "";
      String scans = "";
      int cntr = 0;
      std::vector<FLASHDeconvHelperStructs::Tag> tags;
      std::set<int> sns;
      tagger.getTagsMatchingTo(hit, tags);
      for (const auto& tag : tags)
      {
        if (! tagindices.empty()) tagindices += ";";
        tagindices += std::to_string(tag.getIndex());
        sns.insert(tag.getScan());
        cntr++;
      }

      for (const auto& sn : sns)
      {
        if (! scans.empty()) scans += ";";
        scans += std::to_string(sn);
      }

      fs << tagger.getProteinIndex(hit) << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getSequence() << "\t"
         << hit.getMetaValue("MatchedAA") << "\t" << 100.0 * hit.getCoverage() << "\t" << hit.getScore() << "\t"
         << std::to_string((hit.metaValueExists("qvalue") ? (double)hit.getMetaValue("qvalue") : -1)) << "\t" << cntr << "\t" << tagindices << "\t" << scans << "\n";
      }
    }
  } // namespace OpenMS