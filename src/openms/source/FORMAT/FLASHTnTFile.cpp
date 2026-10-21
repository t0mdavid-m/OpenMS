// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong, Ayesha Feroz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/CHEMISTRY/ProForma.h> // added for ProForma
#include <OpenMS/FORMAT/FLASHTnTFile.h>
#include <exception>
#include <sstream>

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
void FLASHTnTFile::writePrSMHeader(std::fstream& fs, bool chimeric)
{
  fs << "PrSMIndex\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tPrecursorMass\tProteoformMass\tProteoformMassFromComplementaryFragmentIonPairs\tDatabaseSequence\tProteinSequence\tProf"
        "orma\tMatchingFragments\tCoverage(%)\tStartPosition\tEndPosition"
        "\tTagCount\tTagIndices\tModCount\tModMass\tModID\tModAccession\tModStart\tModEnd\tPrecursorQscore\tPrecursorSNR\tScore\tPrSMLevelQvalue\tProteoformLevelQvalue"
        "\tProteoformIndex\tCandidateSource\tFDRClass\tModVariable\tDecoyPrecursorWeight\tProteinLevelQvalue"
        "\tExtensionScore\tRawScore\tRawMatchedFragments\tRawFragmentHypotheses\tRawCoverage(%)"
     << (chimeric ? "\tChimericRound\tMatchedMasses\tRemovedMasses\tResidualNullScore\tResidualScoreRatio" : "") << "\n";
}

/// write header line for Proteoform file
void FLASHTnTFile::writeProHeader(std::fstream& fs, bool chimeric)
{
  fs << "ProteoformIndex\tPrSMIndices\tScan\tRetentionTime\tNumMass\tProteinAccession\tProteinDescription\tPrecursorMass\tProteoformMass\tProteoformMassFromComplementaryFragmentIonPairs\tDatabaseSequence\tProteinSequence\tProf"
        "orma\tMatchingFragments\tCoverage(%)\tStartPosition\tEndPosition"
        "\tTagCount\tTagIndices\tModCount\tModMass\tModID\tModAccession\tModStart\tModEnd\tPrecursorQscore\tPrecursorSNR\tScore\tProteoformLevelQvalue\tCandidateSource\tFDRClass\tModVariable\tProteinLevelQvalue"
        "\tExtensionScore\tRawScore\tRawMatchedFragments\tRawFragmentHypotheses\tRawCoverage(%)"
     << (chimeric ? "\tChimericRound\tMatchedMasses\tRemovedMasses\tResidualNullScore\tResidualScoreRatio" : "") << "\n";
}

/// write the features in regular file output
void FLASHTnTFile::writeProteinHeader(std::fstream& fs)
{
  fs << "Rank\tScan\tProteinAccession\tProteinDescription\tProteinLength\tScore\tMatchedTagCount\tMatchedPositionCount\tTagIndices\n";
}

void FLASHTnTFile::writeProteins(const std::vector<ProteinHit>& hits, std::fstream& fs)
{
  // one line per candidate protein of every spectrum (hundreds of thousands): formatted in parallel blocks, written in order; the meta
  // values by index, because a lookup by name takes the registry's critical section
  auto& registry = MetaInfoInterface::metaRegistry();
  const UInt mi_tag_indices = registry.registerName("TagIndices"), mi_rank = registry.registerName("ProteinRank"), mi_scan = registry.registerName("Scan"),
             mi_matched_aa = registry.registerName("MatchedAA");
  const SignedSize block_size = 2048;
  const SignedSize block_count = ((SignedSize)hits.size() + block_size - 1) / block_size;
  std::vector<std::string> lines(block_count);
  std::vector<std::exception_ptr> errors(block_count);
#pragma omp parallel for schedule(dynamic)
  for (SignedSize b = 0; b < block_count; b++)
  {
    try
    {
      std::stringstream ss; // the default formatting, as that of the file stream
      for (SignedSize i = b * block_size; i < std::min((SignedSize)hits.size(), (b + 1) * block_size); i++)
      {
        const auto& hit = hits[i];
        const std::vector<int> tag_indices = hit.getMetaValue(mi_tag_indices).toIntList();
        ss << (int)hit.getMetaValue(mi_rank) << "\t" << (int)hit.getMetaValue(mi_scan) << "\t" << hit.getAccession() << "\t" << hit.getDescription()
           << "\t" << hit.getSequence().size() << "\t" << hit.getScore() << "\t" << tag_indices.size() << "\t" << (int)hit.getMetaValue(mi_matched_aa) << "\t";
        for (int ti : tag_indices) { ss << ti << ";"; }
        ss << "\n";
      }
      lines[b] = ss.str();
    }
    catch (...)
    {
      errors[b] = std::current_exception();
    }
  }
  // an error is rethrown before anything is written, as by the serial loop
  for (const auto& e : errors)
  {
    if (e) std::rethrow_exception(e);
  }
  for (const auto& l : lines) fs << l;
}

void FLASHTnTFile::writeTags(const FLASHTnTAlgorithm& tnt, double flanking_mass_tol, std::fstream& fs)
{
  const auto& tags = tnt.getTags();
  const auto& proteoforms = tnt.getProteoforms();
  const std::vector<int> no_hits;

  // the lines of the tags with PrSMs first, then those of the others, each in tag order; formatted in parallel by blocks of tags
  const int tag_count = (int)tags.size();
  const int block_size = 4096;
  const int block_count = (tag_count + block_size - 1) / block_size;
  std::vector<std::string> with_hits(block_count), without_hits(block_count);
  std::vector<std::exception_ptr> errors(block_count);
#pragma omp parallel for schedule(dynamic)
  for (int b = 0; b < block_count; b++)
  {
    try
    {
      std::stringstream ss_with, ss_without;
      for (int t = b * block_size; t < std::min(tag_count, (b + 1) * block_size); t++)
      {
        const auto& tag = tags[t];
        const std::vector<int>* matched = tag.getIndex() < 0 ? nullptr : tnt.getProteoformIndicesMatchedBy(tag.getIndex());
        const auto& hit_indices = matched == nullptr ? no_hits : *matched;
        std::stringstream& ss = hit_indices.empty() ? ss_without : ss_with;

        String acc = "";
        String description = "";
        String hitindices = "";
        String positions = "";
        String delta_masses = "";
        for (const int hit_index : hit_indices)
        {
          const auto& hit = proteoforms[hit_index];
          if (! hitindices.empty()) {
            acc += ";";
            description += ";";
            hitindices += ";";
            positions += ";";
            delta_masses += ";";
          }

          acc += hit.getAccession();

          String proteindescription = hit.getDescription();
          if (proteindescription.empty()) { proteindescription = " "; }
          description += proteindescription;
          hitindices += (String)hit.getMetaValue("Index");

          auto pos = std::vector<int>();
          auto masses = std::vector<double>();
          auto pos_in_truncated = std::vector<int>();
          auto masses_in_truncated = std::vector<double>();

          int protein_start_position = hit.getMetaValue("StartPosition");
          int protein_end_position = hit.getMetaValue("EndPosition");
          protein_start_position--;
          String seq = hit.getSequence();
          FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(pos, masses, -1, seq, tag);

          if (protein_end_position >= 0) seq = seq.substr(0, protein_end_position);
          if (protein_start_position >= 0) seq = seq.substr(protein_start_position);
          FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(pos_in_truncated, masses_in_truncated, flanking_mass_tol, seq, tag);
          if (!pos_in_truncated.empty())
          {
            for (int i = 0; i < pos.size(); i++)
            {
              if (pos[i] - (protein_start_position < 0 ? 0 : protein_start_position) != pos_in_truncated[0]) continue;
              positions += std::to_string(pos[i] + 1);
              delta_masses += std::to_string(masses[i]);
              break;
            }
          }
        }

        ss << tag.getIndex() << "\t" << tag.getScan() << "\t" << tag.getRetentionTime() << "\t" << hitindices << "\t" << acc << "\t" << description
           << "\t" << tag.getSequence() << "\t" << std::to_string(tag.getNtermMass()) << "\t" << std::to_string(tag.getCtermMass()) << "\t" << positions
           << "\t" << delta_masses << "\t" << tag.getLength() << "\t" << tag.getScore() << "\t";

        for (const auto& mz : tag.getMzs())
        {
          ss << std::to_string(mz) << ",";
        }
        ss << "\t";
        for (size_t i = 0; i < tag.getMzs().size(); i++) // one score per mass (a gapped tag spans more residues than masses)
        {
          ss << std::to_string(tag.getScore(i)) << ",";
        }
        ss << "\n";
      }
      with_hits[b] = ss_with.str();
      without_hits[b] = ss_without.str();
    }
    catch (...)
    {
      errors[b] = std::current_exception();
    }
  }
  // an error is rethrown before anything is written, as by the serial loop
  for (const auto& e : errors)
  {
    if (e) std::rethrow_exception(e);
  }
  for (const auto& lines : with_hits) fs << lines;
  for (const auto& lines : without_hits) fs << lines;
}

namespace
{
  /// where the candidate protein of a hit came from
  String candidateSource(const OpenMS::ProteinHit& hit)
  {
    if (hit.metaValueExists("PrecursorLadder")) return "precursor ladder";
    if (hit.metaValueExists("TerminalLadder")) return "terminal ladder";
    if (hit.metaValueExists("SlidingLadder")) return "sliding ladder";
    if (hit.metaValueExists("GappedTag")) return "gapped tag";
    return "sequence tag";
  }

  /// 1 for a variable (searched) modification, 0 for a blind shift, per modification
  OpenMS::String modVariable(const OpenMS::ProteinHit& hit)
  {
    OpenMS::String out;
    if (! hit.metaValueExists("ModificationVariable")) return out;
    for (int v : hit.getMetaValue("ModificationVariable").toIntList())
    {
      if (! out.empty()) out += ";";
      out += v != 0 ? "1" : "0";
    }
    return out;
  }

  /// the raw-evidence columns: extension score, raw score, matched fragments, fragment hypotheses, coverage (-1 without a raw file)
  OpenMS::String rawColumns(const OpenMS::ProteinHit& hit)
  {
    OpenMS::String out = "\t";
    out += OpenMS::String(hit.metaValueExists("ExtensionScore") ? (double)hit.getMetaValue("ExtensionScore") : hit.getScore());
    out += "\t" + OpenMS::String(hit.metaValueExists("RawScore") ? (double)hit.getMetaValue("RawScore") : -1.0);
    out += "\t" + OpenMS::String(hit.metaValueExists("RawMatchedFragments") ? (int)hit.getMetaValue("RawMatchedFragments") : -1);
    out += "\t" + OpenMS::String(hit.metaValueExists("RawFragmentHypotheses") ? (int)hit.getMetaValue("RawFragmentHypotheses") : -1);
    out += "\t" + OpenMS::String(hit.metaValueExists("RawCoverage") ? (double)hit.getMetaValue("RawCoverage") : -1.0);
    return out;
  }

  /// the chimeric columns (FLASHTnT -chimeric_rounds): round, matched masses, removed masses, residual null score, score ratio
  OpenMS::String chimericColumns(const OpenMS::ProteinHit& hit)
  {
    OpenMS::String out = "\t";
    out += OpenMS::String(hit.metaValueExists("ChimericRound") ? (int)hit.getMetaValue("ChimericRound") : 0);
    out += "\t" + OpenMS::String(hit.metaValueExists("MatchedMassIndices") ? (int)hit.getMetaValue("MatchedMassIndices").toIntList().size() : -1);
    out += "\t" + OpenMS::String(hit.metaValueExists("RemovedMasses") ? (int)hit.getMetaValue("RemovedMasses") : 0);
    out += "\t" + OpenMS::String(hit.metaValueExists("ResidualNullScore") ? (double)hit.getMetaValue("ResidualNullScore") : -1.0);
    out += "\t" + OpenMS::String(hit.metaValueExists("ResidualScoreRatio") ? (double)hit.getMetaValue("ResidualScoreRatio") : -1.0);
    return out;
  }
} // namespace

String FLASHTnTFile::getProFormaString(const ProteinHit& hit)
{
  const int start = hit.getMetaValue("StartPosition"), end = hit.getMetaValue("EndPosition");
  const String& sequence = hit.getSequence();
  if (start <= 0 || end <= 0 || end < start || end > (int)sequence.size()) return "";
  const int offset = start - 1, length = end - offset;
  ProForma pf;
  try
  {
    pf = ProForma(sequence.substr(offset, length));
  }
  catch (const Exception::ParseError&)
  {
    return ""; // a residue ProForma has no letter for
  }

  struct Entry
  {
    int start, end;
    double mass;
    String names;
    bool variable;
  };
  std::vector<Entry> entries;
  if (hit.metaValueExists("Modifications"))
  {
    const std::vector<double> masses = hit.getMetaValue("Modifications");
    const std::vector<int> starts = hit.getMetaValue("ModificationStarts");
    const std::vector<int> ends = hit.getMetaValue("ModificationEnds");
    const std::vector<String> ids = hit.metaValueExists("ModificationIDs") ? hit.getMetaValue("ModificationIDs").toStringList() : std::vector<String>();
    const std::vector<int> variable = hit.metaValueExists("ModificationVariable") ? (std::vector<int>)hit.getMetaValue("ModificationVariable").toIntList() : std::vector<int>();
    for (Size i = 0; i < masses.size() && i < starts.size() && i < ends.size(); i++)
    {
      const int s = std::max(0, std::min(length - 1, std::min(starts[i], ends[i]) - offset));
      const int e = std::max(0, std::min(length - 1, std::max(starts[i], ends[i]) - offset));
      entries.push_back({s, e, masses[i], i < ids.size() ? ids[i] : String(), i < variable.size() && variable[i] != 0});
    }
  }
  std::stable_sort(entries.begin(), entries.end(), [](const Entry& l, const Entry& r) { return l.start < r.start; });
  // ProForma ranges are the same or disjoint: two ranges that overlap otherwise become one, with their summed mass
  std::vector<Entry> merged;
  for (const auto& entry : entries)
  {
    bool combined = false;
    for (auto& m : merged)
    {
      if (entry.end == entry.start || m.end == m.start || (m.start == entry.start && m.end == entry.end)) continue;
      if (m.start <= entry.end && entry.start <= m.end)
      {
        m.start = std::min(m.start, entry.start);
        m.end = std::max(m.end, entry.end);
        m.mass += entry.mass;
        m.names = "";
        m.variable = false;
        combined = true;
        break;
      }
    }
    if (! combined) merged.push_back(entry);
  }
  for (const auto& m : merged)
  {
    // the names of a variable modification ("A," or "A+B,"; alternatives "A,B," leave the identity open)
    std::vector<std::pair<String, double>> named;
    String names = m.names;
    while (names.hasSuffix(",")) names.pop_back();
    if (m.variable && ! names.empty() && ! names.has(','))
    {
      std::vector<String> parts;
      names.split('+', parts);
      double sum = 0;
      for (const auto& part : parts)
      {
        try
        {
          const double mass = ProForma::getTagMass(part, pf.getSequence()[m.start]);
          named.emplace_back(part, mass);
          sum += mass;
        }
        catch (const Exception::ParseError&)
        {
          named.clear();
          break;
        }
      }
      if (std::abs(sum - m.mass) > 0.05) named.clear();
    }
    try
    {
      if (named.empty()) pf.addModification(m.start, m.end, m.mass);
      for (const auto& [name, mass] : named) pf.addModification(m.start, m.end, mass, name);
    }
    catch (const Exception::InvalidValue&)
    {
      return "";
    }
  }
  return pf.toString();
}


void OpenMS::FLASHTnTFile::writePrSMs(const std::vector<ProteinHit>& hits, std::fstream& fs, bool chimeric)
{
  std::stringstream ss;
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
      if (i > 0)
      {
        modmasses += ";";
        modstarts += ";";
        modends += ";";
        modids += ";";
        modaccs += ";";
      }

      modmasses += std::to_string(mod_masses[i]);
      modstarts += std::to_string(mod_starts[i] + 1);
      modends += std::to_string(mod_ends[i] + 1);
      modids += mod_ids[i];
      modaccs += mod_accs[i];
    }

    int start = hit.getMetaValue("StartPosition");
    int end = hit.getMetaValue("EndPosition");

    int start_in_seq = start < 0 ? 0 : (start - 1);
    int end_in_seq = end < 0 ? hit.getSequence().size() : end;

    // Use ProForma for sequence generation
    String proformaStr = getProFormaString(hit);
    ss << hit.getMetaValue("Index") << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getMetaValue("RT") << "\t" << hit.getMetaValue("NumMass")
       << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getMetaValue("GivenMass")  << "\t" << hit.getMetaValue("Mass") << "\t" << ((int)hit.getMetaValue("ProteoformMassByFragmentMass") > 0) << "\t" <<  hit.getSequence() << "\t"
       << hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << hit.getMetaValue("MatchedAA") << "\t"
       << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << cntr << "\t" << tagindices << "\t" << mod_masses.size() << "\t"
       << modmasses << "\t" << modids << "\t" << modaccs << "\t" << modstarts << "\t" << modends << "\t" << hit.getMetaValue("PrecursorScore")  << "\t" << hit.getMetaValue("PrecursorSNR")  << "\t" << hit.getScore() << "\t"
       << std::to_string((hit.metaValueExists("qvalue") ? (double)hit.getMetaValue("qvalue") : -1)) << "\t"
       << std::to_string((hit.metaValueExists("proqvalue") ? (double)hit.getMetaValue("proqvalue") : -1)) << "\t"
       << (hit.metaValueExists("ProteoformIndex") ? (int)hit.getMetaValue("ProteoformIndex") : -1) << "\t" << candidateSource(hit) << "\t"
       << (hit.metaValueExists("FDRClass") ? hit.getMetaValue("FDRClass").toString() : String()) << "\t" << modVariable(hit) << "\t"
       << (hit.metaValueExists("DecoyPrecursor") ? (double)hit.getMetaValue("DecoyPrecursor") : 0.0) << "\t"
       << std::to_string(hit.metaValueExists("ProteinQvalue") ? (double)hit.getMetaValue("ProteinQvalue") : -1) << rawColumns(hit)
       << (chimeric ? chimericColumns(hit) : String()) << "\n";
  }
  fs << ss.str();
}

void OpenMS::FLASHTnTFile::writeProteoforms(const std::vector<ProteinHit>& hits, std::fstream& fs, double pro_fdr, bool chimeric)
{
  std::stringstream ss;
  int pro_index = 0;
  for (const auto& hit : hits)
  {
    if (! hit.metaValueExists("Index")) continue;
    if ((int)hit.getMetaValue("Representative") == 0) continue;
    if (hit.metaValueExists("proqvalue") && (double)hit.getMetaValue("proqvalue") > pro_fdr) continue;
    String prsmindices = "", tagindices = "", modmasses = "", modstarts = "", modends = "", modids = "", modaccs = "";

    int cntr = 0;
    //std::vector<FLASHHelperClasses::Tag> tags;
    std::vector<int> tag_indices;
    if (hit.metaValueExists("TagIndices")) tag_indices = (std::vector<int>)hit.getMetaValue("TagIndices").toIntList();
    for (const int& index : tag_indices)
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
      if (i > 0)
      {
        modmasses += ";";
        modstarts += ";";
        modends += ";";
        modids += ";";
        modaccs += ";";
      }

      modmasses += std::to_string(mod_masses[i]);
      modstarts += std::to_string(mod_starts[i] + 1);
      modends += std::to_string(mod_ends[i] + 1);
      modids += mod_ids[i];
      modaccs += mod_accs[i];
    }
    if (hit.metaValueExists("PrSMIndices"))
    {
      auto prsm_indices = hit.getMetaValue("PrSMIndices").toIntList();
      for (const int& index : prsm_indices)
      {
        if (! prsmindices.empty()) prsmindices += ";";
        prsmindices += std::to_string(index);
      }
    }
    int start = hit.getMetaValue("StartPosition");
    int end = hit.getMetaValue("EndPosition");

    int start_in_seq = start < 0 ? 0 : (start - 1);
    int end_in_seq = end < 0 ? hit.getSequence().size() : end;

    // Use ProForma to generate ProForma string
    String proformaStr = getProFormaString(hit);

    const int this_index = pro_index++;
    ss << (hit.metaValueExists("ProteoformIndex") ? (int)hit.getMetaValue("ProteoformIndex") : this_index) << "\t" << prsmindices << "\t" << hit.getMetaValue("Scan") << "\t" << hit.getMetaValue("RT") << "\t" << hit.getMetaValue("NumMass")
       << "\t" << hit.getAccession() << "\t" << hit.getDescription() << "\t" << hit.getMetaValue("GivenMass") << "\t" << hit.getMetaValue("Mass") << "\t" << ((int)hit.getMetaValue("ProteoformMassByFragmentMass") > 0) << "\t" << hit.getSequence() << "\t"
       << hit.getSequence().substr(start_in_seq, end_in_seq - start_in_seq) << "\t" << proformaStr << "\t" << hit.getMetaValue("MatchedAA") << "\t"
       << 100.0 * hit.getCoverage() << "\t" << start << "\t" << end << "\t" << cntr << "\t" << tagindices << "\t" << mod_masses.size() << "\t"
       << modmasses << "\t" << modids << "\t" << modaccs << "\t" << modstarts << "\t" << modends << "\t" << hit.getMetaValue("PrecursorScore") << "\t" << hit.getMetaValue("PrecursorSNR") << "\t" << hit.getScore() << "\t"
       << std::to_string((hit.metaValueExists("proqvalue") ? (double)hit.getMetaValue("proqvalue") : -1)) << "\t" << candidateSource(hit) << "\t"
       << (hit.metaValueExists("FDRClass") ? hit.getMetaValue("FDRClass").toString() : String()) << "\t" << modVariable(hit) << "\t"
       << std::to_string(hit.metaValueExists("ProteinQvalue") ? (double)hit.getMetaValue("ProteinQvalue") : -1) << rawColumns(hit)
       << (chimeric ? chimericColumns(hit) : String()) << "\n";
  }
  fs << ss.str();
}

} // namespace OpenMS
// namespace OpenMS