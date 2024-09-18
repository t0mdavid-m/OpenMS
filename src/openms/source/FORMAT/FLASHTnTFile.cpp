// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong, Ayesha Feroz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/FORMAT/FLASHTnTFile.h>
#include <map>

namespace OpenMS
{
/**
  @brief FLASHTnT output *.tsv file format
   @ingroup FileIO
**/
    void FLASHTnTFile::writeTagHeader(std::fstream& fs)
    {
        fs << "TagIndex\tScan\tRetentionTime\tProteoformIndex\tProteinAccession\tProteinDescription\tTagSequence\tNmass\tCmass\tStartPosition\tDeltaMass\tLength\tDeNovoScore\tMasses\tMassScores\n";
    }

/// write header line for PrSM file
    void FLASHTnTFile::writePrSMHeader(std::fstream& fs)
    {
        fs << "ProteoformIndex\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tProteoformMass\tDatabaseSequence\tProteinSequence\tProforma\tMatchedAminoAcidCount\tCoverage(%)\tStartPosition\tEndPosition"
              "\tTagCount\tTagIndices\tModCount\tModMass\tModID\tModAccession\tModStart\tModEnd\tScore\tPrSMLevelQvalue\tProteoformLevelQvalue\n";
    }

/// write header line for Proteoform file
    void FLASHTnTFile::writeProHeader(std::fstream& fs)
    {
        fs << "ProteoformIndex\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tProteoformMass\tDatabaseSequence\tProteinSequence\tProforma\tMatchedAminoAcidCount\tCoverage(%)\tStartPosition\tEndPosition"
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
                if (c == 1 && !hits.empty()) continue;

                String acc = "";
                String description = "";
                String hitindices = "";
                String positions = "";
                String delta_masses = "";
                for (const auto& hit : hits)
                {
                    if (!acc.empty()) acc += ";";
                    if (!description.empty()) description += ";";
                    if (!hitindices.empty()) hitindices += ";";
                    if (!positions.empty()) positions += ";";
                    if (!delta_masses.empty()) delta_masses += ";";

                    acc += hit.getAccession();

                    String proteindescription = hit.getDescription();
                    if (proteindescription.empty())
                    {
                        proteindescription = " ";
                    }
                    description += proteindescription;
                    hitindices += (String)hit.getMetaValue("Index");

                    auto pos = std::vector<int>();
                    auto masses = std::vector<double>();
                    FLASHTaggerAlgorithm::getMatchedPositionsAndFlankingMassDiffs(pos, masses, flanking_mass_tol, hit, tag);
                    if (pos.size() != 0)
                    {
                        positions += std::to_string(pos[0]);
                    }
                    if (masses.size() != 0)
                    {
                        delta_masses += std::to_string(masses[0]);
                    }
                }

                fs << tag.getIndex() << "\t" << tag.getScan() << "\t" << tag.getRetentionTime() << "\t" << hitindices << "\t" << acc << "\t" << description
                   << "\t" << tag.getSequence() << "\t" << std::to_string(tag.getNtermMass()) << "\t" << std::to_string(tag.getCtermMass()) << "\t" << positions
                   << "\t" << delta_masses << "\t" << tag.getLength() << "\t" << tag.getScore() << "\t";

                for (const auto& mz : tag.getMzs())
                {
                    fs << std::to_string(mz) << ",";
                }
                fs << "\t";
                for (auto i = 0; i <= tag.getLength(); i++)
                {
                    fs << std::to_string(tag.getScore(i)) << ",";
                }
                fs << "\n";
            }
        }
    }

// Helper function to generate ProForma string
    String generateProFormaString(const String& sequence, const std::vector<double>& mod_masses, const std::vector<int>& mod_starts, const std::vector<int>& mod_ends, const std::vector<String>& mod_ids, const std::vector<String>& mod_accs)
    {
        String proformaStr = sequence;
        std::map<int, String> modifications_map;

        for (size_t i = 0; i < mod_masses.size(); ++i)
        {
            String mod_str = "[";

            if (!mod_ids[i].empty())
            {
                mod_str += mod_ids[i];
            }
            else if (mod_masses[i] != 0.0)
            {
                mod_str += (mod_masses[i] > 0 ? "+" : "") + std::to_string(mod_masses[i]);
            }

            if (!mod_accs[i].empty())
            {
                mod_str += ":" + mod_accs[i];
            }

            mod_str += "]";

            // Apply the modification at the correct position
            int position = mod_starts[i] + 1; // ProForma is 1-based
            modifications_map[position] = mod_str;
        }

        // Now, integrate modifications into the sequence
        String proforma_sequence;
        for (size_t i = 0; i < sequence.size(); ++i)
        {
            proforma_sequence += sequence[i];
            if (modifications_map.count(i + 1) > 0)
            {
                proforma_sequence += modifications_map[i + 1];
            }
        }

        return proforma_sequence;
    }

    void FLASHTnTFile::writePrSMs(const std::vector<ProteinHit>& hits, std::fstream& fs)
    {
        for (const auto& hit : hits)
        {
            if (!hit.metaValueExists("Index")) continue;

            // Extract modifications
            std::vector<double> mod_masses = hit.getMetaValue("Modifications");
            std::vector<int> mod_starts = hit.getMetaValue("ModificationStarts");
            std::vector<int> mod_ends = hit.getMetaValue("ModificationEnds");
            std::vector<String> mod_ids = hit.getMetaValue("ModificationIDs");
            std::vector<String> mod_accs = hit.getMetaValue("ModificationACCs");

            int start = hit.getMetaValue("StartPosition");
            int end = hit.getMetaValue("EndPosition");

            int start_in_seq = start < 0 ? 0 : (start - 1);
            int end_in_seq = end < 0 ? hit.getSequence().size() : end;

            // Generate ProForma string
            String proformaStr = generateProFormaString(hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq), mod_masses, mod_starts, mod_ends, mod_ids, mod_accs);

            fs << hit.getMetaValue("Index") << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getMetaValue("RT") << "\t" << hit.getMetaValue("NumMass") << "\t" << hit.getAccession() << "\t"
               << hit.getDescription() << "\t" << hit.getMetaValue("Mass") << "\t" << hit.getSequence() << "\t"
               << hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << hit.getMetaValue("MatchedAA") << "\t"
               << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << mod_masses.size() << "\t" << hit.getScore() << "\t"
               << std::to_string((hit.metaValueExists("qvalue") ? (double)hit.getMetaValue("qvalue") : -1)) << "\n";
        }
    }

    void FLASHTnTFile::writeProteoforms(const std::vector<ProteinHit>& hits, std::fstream& fs, double pro_fdr)
    {
        for (const auto& hit : hits)
        {
            if (!hit.metaValueExists("Index")) continue;
            if (!hit.metaValueExists("Representative")) continue;
            if (hit.metaValueExists("proqvalue") && (double) hit.getMetaValue("proqvalue") > pro_fdr) continue;

            // Extract modifications
            std::vector<double> mod_masses = hit.getMetaValue("Modifications");
            std::vector<int> mod_starts = hit.getMetaValue("ModificationStarts");
            std::vector<int> mod_ends = hit.getMetaValue("ModificationEnds");
            std::vector<String> mod_ids = hit.getMetaValue("ModificationIDs");
            std::vector<String> mod_accs = hit.getMetaValue("ModificationACCs");

            int start = hit.getMetaValue("StartPosition");
            int end = hit.getMetaValue("EndPosition");

            int start_in_seq = start < 0 ? 0 : (start - 1);
            int end_in_seq = end < 0 ? hit.getSequence().size() : end;

            // Generate ProForma string
            String proformaStr = generateProFormaString(hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq), mod_masses, mod_starts, mod_ends, mod_ids, mod_accs);

            fs << hit.getMetaValue("Index") << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getMetaValue("RT") << "\t" << hit.getMetaValue("NumMass") << "\t" << hit.getAccession() << "\t"
               << hit.getDescription() << "\t" << hit.getMetaValue("Mass") << "\t" << hit.getSequence() << "\t"
               << hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << hit.getMetaValue("MatchedAA") << "\t"
               << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << mod_masses.size() << "\t" << hit.getScore() << "\t"
               << std::to_string((hit.metaValueExists("proqvalue") ? (double)hit.getMetaValue("proqvalue") : -1)) << "\n";
        }
    }

}
 // namespace OpenMS