// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <fstream>
#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_in_parallel() 0
#endif
#include <queue>
#include <utility>

namespace OpenMS
  {
  inline const int multi_ion_score = 1;
  inline const bool debug = false;
  inline const double i2f_mass = Residue::getInternalToFull().getMonoWeight();
  // Base penalty for a terminal truncation, equivalent to the base blind-modification (mass-shift)
  // penalty applied in findSubPathsBetweenTagEndPoints. Equals 19.
  inline const int truncation_penalty = 1 + 2 * (multi_ion_score + FLASHTaggerAlgorithm::max_node_score);
  // Minimum number of unshifted suffix fragment masses from the protein C-terminus against which a precursor-mass path ending before
  // the C-terminus is not accepted: terminal ladders of k >= 6 had no chance hits in the own-candidate entrapment study
  // (final_all/sliding/own_cand).
  inline const int min_c_term_ladder = 6;
  // Penalty per variable modification: a known residue-specific shift costs far less than a blind modification, but more than
  // nothing, so that an unmodified path wins a tie.
  inline const int var_mod_penalty = 1 + multi_ion_score;
  FLASHExtenderAlgorithm::FLASHExtenderAlgorithm(): DefaultParamHandler("FLASHExtenderAlgorithm"), ProgressLogger()
  {
    setDefaultParams_();
  }

  FLASHExtenderAlgorithm& FLASHExtenderAlgorithm::operator=(const FLASHExtenderAlgorithm& rhs)
  {
    if (this == &rhs) return *this;

    DefaultParamHandler::operator=(rhs);
    return *this;
  }

  void FLASHExtenderAlgorithm::setDefaultParams_()
  {
    defaults_.setValue("max_mod_mass", 500.0, "Maximum mass shift allowed for modifications.");
    defaults_.setValue("max_mod_count", 2, "Maximum number of blind modifications.");

    defaults_.setValue("var_mod", std::vector<std::string> {}, "Variable modifications, by name as in Unimod/PSI-MOD (e.g. 'Oxidation (M)', "
                       "'Acetyl (N-term)', 'Phospho (S)'). Searched with their own count, independently of blind modifications; both can occur "
                       "in one proteoform.");
    defaults_.setValue("var_mod_file", "", "Variable modifications in TopPIC's format, one per line: name,mass,residues,position,unimod "
                       "(residues: one-letter codes or *; position: any, N-term or C-term; # starts a comment). Added to var_mod.");
    defaults_.setValue("max_var_mod_count", 3, "Maximum number of variable modifications per proteoform (used only when variable modifications "
                       "are given).");
    defaults_.setMinInt("max_var_mod_count", 0);
    defaults_.setValue("var_mod_penalty", std::vector<int> {},
                       "Score penalty of each variable modification, in the order given (var_mod, then var_mod_file); missing or lower: 2. "
                       "FLASHTnT raises it for modifications its decoy PrSMs use about as often as its target PrSMs.");
    defaults_.addTag("var_mod_penalty", "advanced");

    defaults_.setValue("cleavage_weight", 5.0, "Weight of the backbone cleavage propensity: a matched fragment mass adds round(weight x (w_before + "
                       "w_after)) to the path score (at most +-4), where w is log2(cleavage rate / mean) of the residues before and after the bond "
                       "for its series (b/y or c/z). 0 switches it off.");
    defaults_.setMinFloat("cleavage_weight", 0.0);
    defaults_.setValue("cleavage_table", "", "Cleavage propensities replacing the built-in ones, one per line: series (by or cz), before or after, "
                       "residue, log2(rate / mean). Built-in: b/y from HCD E. coli lysates, c/z from ECD native proteins.");

    defaults_.setValue("ion_type", std::vector<std::string> {"b", "y"}, "Specifies ion types to consider");
    defaults_.setValidStrings("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});
    defaults_.setValue("skip_precursor_inference", "false", "Set to skip precursor mass inference with complimentary fragment ion pairs. Precursor masses inferred from MSn-1 will be used in MSn");
    defaults_.setValidStrings("skip_precursor_inference", {"true", "false"});
    defaults_.addTag("skip_precursor_inference", "advanced");

    defaultsToParam_();
  }

  void FLASHExtenderAlgorithm::updateMembers_()
  {
    max_blind_mod_cntr_ = param_.getValue("max_mod_count");
    max_mod_mass_ = param_.getValue("max_mod_mass");
    skip_precursor_inference_ = param_.getValue("skip_precursor_inference") == "true";
    loadVariableModifications_();
    cleavage_weight_ = param_.getValue("cleavage_weight");
    loadCleavageTables_();

    // precomputed vertex-packing spans (invariant during run(); used by the get*Vertex helpers on the hot path)
    var_span_ = (Size)max_var_mod_cntr_ + 1;
    score_span_ = (Size)(max_path_score_ - min_path_score_) + 1;
    blind_span_ = (Size)max_blind_mod_cntr_ + 1;
    var_score_span_ = var_span_ * score_span_;
    var_score_blind_span_ = var_score_span_ * blind_span_;
  }

  void FLASHExtenderAlgorithm::loadCleavageTables_()
  {
    // log2(cleavage rate / mean) of the residue before / after a bond, A C D E F G H I K L M N P Q R S T V W Y.
    // b/y: Eclip_1476/1479/1480 (HCD), 4592 unmodified PrSMs at 1% FDR, 618,683 bonds; the three runs agree (r >= 0.986).
    // c/z: ECD native proteins (NIST mAb, Infliximab, AQPZ, DLDH, OMPA/OMPF/G3P1), 49 PrSMs, 18,005 bonds, shrunk with 600 pseudo-bonds.
    static const char aa[] = "ACDEFGHIKLMNPQRSTVWY";
    static const double by_before[] = {+0.04, -0.90, +0.63, +0.51, -0.25, -0.79, +0.14, +0.33, -0.18, +0.07, -0.56, -0.05, -1.13, +0.13, -1.61, -0.12, -0.07, +0.48, -0.08, -0.23};
    static const double by_after[] = {+0.06, -0.26, -0.07, -0.07, -0.08, -0.05, -0.10, +0.22, -0.42, +0.02, -0.49, +0.01, +0.80, +0.11, -2.09, +0.36, +0.14, +0.25, -0.15, -0.18};
    static const double cz_before[] = {-0.24, +0.22, +0.19, -0.28, +0.10, +0.31, +0.19, -0.18, -0.19, -0.27, -0.29, -0.22, -0.04, +0.14, -0.39, +0.53, +0.14, +0.07, -0.20, -0.11};
    static const double cz_after[] = {-0.24, -0.03, -0.06, -0.30, +0.30, +0.18, +0.09, -0.12, +0.13, -0.24, -0.23, -0.02, -0.45, +0.14, -0.30, +0.49, +0.29, -0.04, -0.05, +0.04};
    for (int s = 0; s < 2; s++)
    {
      for (int r = 0; r < 26; r++) cleavage_before_[s][r] = cleavage_after_[s][r] = 0;
    }
    for (int i = 0; i < 20; i++)
    {
      const int r = aa[i] - 'A';
      cleavage_before_[0][r] = by_before[i];
      cleavage_after_[0][r] = by_after[i];
      cleavage_before_[1][r] = cz_before[i];
      cleavage_after_[1][r] = cz_after[i];
    }
    const String file = param_.getValue("cleavage_table").toString();
    if (file.empty()) return;
    std::ifstream in(file);
    if (! in.good()) throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file);
    String line;
    while (std::getline(in, line))
    {
      line.trim();
      if (line.empty() || line.hasPrefix("#")) continue;
      std::vector<String> f;
      line.split('\t', f);
      if (f.size() < 4 || f[2].size() != 1 || f[2][0] < 'A' || f[2][0] > 'Z' || (f[0] != "by" && f[0] != "cz") || (f[1] != "before" && f[1] != "after"))
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "expected: by|cz <tab> before|after <tab> residue <tab> log2 ratio");
      const int series = f[0] == "by" ? 0 : 1;
      (f[1] == "before" ? cleavage_before_ : cleavage_after_)[series][f[2][0] - 'A'] = f[3].toDouble();
    }
  }

  void FLASHExtenderAlgorithm::loadVariableModifications_()
  {
    var_mod_list_.clear();
    var_mod_combinations_.clear();
    max_var_mod_mass_ = 0;

    const auto inst = ModificationsDB::getInstance();
    for (const auto& name : param_.getValue("var_mod").toStringVector())
    {
      if (String(name).trim().empty()) continue;
      const ResidueModification* mod = inst->getModification(name);
      VariableModification vm;
      vm.name = mod->getId();
      vm.mass = mod->getDiffMonoMass();
      const char origin = mod->getOrigin();
      vm.residues = (origin == 'X' || origin == '\0' || origin == '.') ? String("") : String(1, origin);
      const auto term = mod->getTermSpecificity();
      vm.position = (term == ResidueModification::N_TERM || term == ResidueModification::PROTEIN_N_TERM) ? 1
                    : (term == ResidueModification::C_TERM || term == ResidueModification::PROTEIN_C_TERM) ? 2 : 0;
      vm.unimod = mod->getUniModRecordId();
      var_mod_list_.push_back(vm);
    }

    const String file = param_.getValue("var_mod_file").toString();
    if (! file.empty())
    {
      std::ifstream in(file);
      if (! in.good()) throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file);
      String line;
      while (std::getline(in, line))
      {
        line.trim();
        if (line.empty() || line.hasPrefix("#")) continue;
        std::vector<String> f;
        line.split(',', f);
        if (f.size() < 4) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, line, "expected name,mass,residues,position[,unimod]");
        VariableModification vm;
        vm.name = f[0].trim();
        vm.mass = f[1].toDouble();
        vm.residues = f[2].trim() == "*" ? String("") : f[2].trim();
        if (vm.residues.size() == 2 && vm.residues[0] == vm.residues[1]) // e.g. CC: a disulfide needs two cysteines
        {
          vm.residues = vm.residues.substr(0, 1);
          vm.sites = 2;
        }
        const String pos = String(f[3]).trim().toLower();
        vm.position = pos.hasSubstring("n-term") ? 1 : pos.hasSubstring("c-term") ? 2 : 0;
        if (f.size() > 4 && ! f[4].trim().empty()) vm.unimod = f[4].toInt();
        var_mod_list_.push_back(vm);
      }
    }

    // a penalty per modification; never below the default, so a raised penalty only removes paths
    var_mod_penalties_.assign(var_mod_list_.size(), var_mod_penalty);
    const auto given_penalties = param_.getValue("var_mod_penalty").toIntVector();
    for (Size i = 0; i < std::min(given_penalties.size(), var_mod_penalties_.size()); i++)
    {
      var_mod_penalties_[i] = std::max(var_mod_penalty, given_penalties[i]);
    }

    max_var_mod_cntr_ = var_mod_list_.empty() ? 0 : (int)param_.getValue("max_var_mod_count");
    if (max_var_mod_cntr_ == 0) return;

    for (const auto& vm : var_mod_list_) max_var_mod_mass_ = std::max(max_var_mod_mass_, std::abs(vm.mass));
    // every multiset of up to max_var_mod_cntr_ modifications that can sit on one proteoform (at most one N- and one C-terminal
    // modification), sorted by mass
    std::vector<int> idx;
    var_mod_room_.assign(max_var_mod_cntr_ + 1, 0.0);
    std::function<void(int, double, int, int)> add = [&](int from, double mass, int n_terms, int c_terms) {
      if (! idx.empty())
      {
        var_mod_combinations_.push_back({mass, idx});
        for (int k = (int)idx.size(); k <= max_var_mod_cntr_; k++) var_mod_room_[k] = std::max(var_mod_room_[k], std::abs(mass));
      }
      if ((int)idx.size() == max_var_mod_cntr_) return;
      for (int i = from; i < (int)var_mod_list_.size(); i++)
      {
        const int n = n_terms + (var_mod_list_[i].position == 1), c = c_terms + (var_mod_list_[i].position == 2);
        if (n > 1 || c > 1) continue;
        idx.push_back(i);
        add(i, mass + var_mod_list_[i].mass, n, c);
        idx.pop_back();
      }
    };
    add(0, 0, 0, 0);
    std::sort(var_mod_combinations_.begin(), var_mod_combinations_.end(),
              [](const VariableModificationCombination& a, const VariableModificationCombination& b) { return a.mass < b.mass || (a.mass == b.mass && a.mods < b.mods); });
  }

  bool FLASHExtenderAlgorithm::combinationFits_(const VariableModificationCombination& combination, const String& seq, int from, int to, int n_site,
                                                int c_site) const
  {
    // called in the innermost loop of the extension: no allocation
    to = std::min(to, (int)seq.size());
    auto fits_at = [&](const VariableModification& vm, int p) {
      if (vm.position == 1 && p != n_site) return false; // the proteoform's first residue
      if (vm.position == 2 && p != c_site) return false; // the proteoform's last residue
      return vm.residues.empty() || vm.residues.find(seq[p]) != std::string::npos;
    };
    const auto& mods = combination.mods;
    if (mods.size() == 1 && var_mod_list_[mods[0]].sites == 1)
    {
      const auto& vm = var_mod_list_[mods[0]];
      for (int p = from; p < to; p++)
        if (fits_at(vm, p)) return true;
      return false;
    }
    // distinct residues for all modification sites (a two-site modification takes two): depth-first
    int site_mod[16];
    int k_max = 0;
    for (int m : mods)
    {
      for (int t = 0; t < var_mod_list_[m].sites; t++)
      {
        if (k_max == 16) return false;
        site_mod[k_max++] = m;
      }
    }
    int site[16];
    int k = 0;
    site[0] = from - 1;
    while (k >= 0)
    {
      int p = site[k] + 1;
      for (; p < to; p++)
      {
        if (! fits_at(var_mod_list_[site_mod[k]], p)) continue;
        bool taken = false;
        for (int q = 0; q < k && ! taken; q++) taken = site[q] == p;
        if (! taken) break;
      }
      if (p >= to) { k--; continue; } // backtrack
      site[k] = p;
      if (++k == k_max) return true;
      site[k] = from - 1;
    }
    return false;
  }

  String FLASHExtenderAlgorithm::variableModificationNames_(double mass, double tol, const String& seq, int from, int to, int n_site, int c_site,
                                                            String& accessions, int& site) const
  {
    String names;
    accessions = "";
    site = -1;
    int fitting = 0;
    const VariableModificationCombination* last = nullptr;
    auto iter = std::lower_bound(var_mod_combinations_.begin(), var_mod_combinations_.end(), mass - tol,
                                 [](const VariableModificationCombination& c, double m) { return c.mass < m; });
    for (; iter != var_mod_combinations_.end() && iter->mass <= mass + tol; ++iter)
    {
      if (! combinationFits_(*iter, seq, from, to, n_site, c_site)) continue;
      fitting++;
      last = &*iter;
      for (Size k = 0; k < iter->mods.size(); k++)
      {
        const auto& vm = var_mod_list_[iter->mods[k]];
        names += (k ? "+" : "") + vm.name;
        accessions += (k ? String("+") : String("")) + (vm.unimod >= 0 ? String(vm.unimod) : String("-"));
      }
      names += ",";
      accessions += ",";
    }
    if (fitting == 1 && last->mods.size() == 1 && var_mod_list_[last->mods[0]].sites == 1)
    {
      const auto& vm = var_mod_list_[last->mods[0]];
      int count = 0;
      for (int p = from; p < std::min(to, (int)seq.size()); p++)
      {
        if (vm.position == 1 && p != n_site) continue;
        if (vm.position == 2 && p != c_site) continue;
        if (! vm.residues.empty() && vm.residues.find(seq[p]) == std::string::npos) continue;
        site = p;
        count++;
      }
      if (count != 1) site = -1;
    }
    return names;
  }

  int FLASHExtenderAlgorithm::shiftRangeEnd_(const MSSpectrum& node_spec, const MSSpectrum& tol_spec, const std::vector<double>& pro_masses, int pre_node,
                                             int pre_pro, int node, int pro, double mass_shift) const
  {
    // A mass between the two nodes that matches protein index pi with the whole shift puts the shift on residues before pi; the smallest
    // such pi ends the range. Called on DAG edges (the innermost loop of the extension): per mass in between, the candidate indices are
    // found by binary search in the protein masses (non-decreasing), and each is tested exactly as |(mass_shift + protein mass) - node
    // mass| <= tolerance. The search bounds are 1e-6 Da wider than the tolerance, far above rounding, so no match is skipped.
    int end = pro;
    for (int ni = pre_node + 1; ni < node && end > pre_pro + 1; ni++)
    {
      const double nm = node_spec[ni].getMZ();
      const double margin = tol_spec[ni].getIntensity();
      const auto last = pro_masses.begin() + end; // only indices below the current end can end the range earlier
      auto it = std::lower_bound(pro_masses.begin() + (pre_pro + 1), last, nm - margin - mass_shift - 1e-6);
      for (; it != last && *it <= nm + margin - mass_shift + 1e-6; ++it)
      {
        const double pm = mass_shift + *it;
        if (std::abs(pm - nm) > margin) continue;
        end = (int)(it - pro_masses.begin());
        break;
      }
    }
    return end;
  }

  inline Size FLASHExtenderAlgorithm::getVertex_(int node_index, int pro_index, int score, int num_blind_mod, int num_var_mod, Size pro_mass_size) const
  {
    return (((node_index * pro_mass_size + pro_index) * blind_span_ + num_blind_mod) * score_span_
            + (std::min(max_path_score_, std::max(min_path_score_, score)) - min_path_score_))
             * var_span_
           + num_var_mod;
  }

  inline int FLASHExtenderAlgorithm::getNodeIndex_(Size vertex, Size pro_mass_size) const
  {
    return vertex / var_score_blind_span_ / pro_mass_size;
  }

  inline int FLASHExtenderAlgorithm::getProIndex_(Size vertex, Size pro_mass_size) const
  {
    return (vertex / var_score_blind_span_) % pro_mass_size;
  }

  inline int FLASHExtenderAlgorithm::getScore_(Size vertex) const
  {
    return (vertex / var_span_) % score_span_ + min_path_score_;
  }

  inline int FLASHExtenderAlgorithm::getBlindModNumber_(Size vertex) const
  {
    return (vertex / var_score_span_) % blind_span_;
  }

  inline int FLASHExtenderAlgorithm::getVarModNumber_(Size vertex) const
  {
    return vertex % var_span_;
  }

  // take the hits. Just calculate the mass of truncated protein. Then add modification masses if they are disjoint. If they overlap and the same mass,
  // we have a single one. If they are overlapping but different, we add all of them. the max mod count is also adjusted.
  bool FLASHExtenderAlgorithm::hasComplementaryPairSupport_(const std::vector<double>& sorted_masses, double precursor_mass) const
  {
    if (precursor_mass <= 0 || sorted_masses.size() < 2) return false;
    // a prefix and a suffix fragment of one precursor sum to the precursor mass plus their ion offsets minus water
    std::vector<double> offsets;
    for (double p : prefix_shifts_)
      for (double q : suffix_shifts_)
        offsets.push_back(p + q - i2f_mass);
    if (offsets.empty()) return false;
    auto count_pairs = [&](double mass) {
      int count = 0;
      const double margin = tol_ * mass;
      for (double offset : offsets)
      {
        const double target = mass + offset;
        for (Size i = 0; i < sorted_masses.size() && sorted_masses[i] * 2 < target + margin; i++)
        {
          auto lo = std::lower_bound(sorted_masses.begin() + i + 1, sorted_masses.end(), target - sorted_masses[i] - margin);
          for (; lo != sorted_masses.end() && *lo <= target - sorted_masses[i] + margin; ++lo) count++;
        }
      }
      return count;
    };
    const int observed = count_pairs(precursor_mass);
    if (observed == 0) return false;
    // Chance: the same count at masses moved by up to +-300 Da, avoiding shifts that real neighbours of the precursor produce (isotopes,
    // water / ammonia loss, oxidation, sodium). In a spectrum crowded with masses (native complexes, 200-400 masses) the moved masses
    // collect as many pairs as the inferred one; in a sparse one (lysate MS2) they collect almost none.
    static const std::vector<double> neighbours {-4 * Constants::C13C12_MASSDIFF_U, -3 * Constants::C13C12_MASSDIFF_U, -2 * Constants::C13C12_MASSDIFF_U,
                                                 -Constants::C13C12_MASSDIFF_U,     0.0,                                Constants::C13C12_MASSDIFF_U,
                                                 2 * Constants::C13C12_MASSDIFF_U,  3 * Constants::C13C12_MASSDIFF_U,  4 * Constants::C13C12_MASSDIFF_U,
                                                 -18.010565,                        -17.026549,                         15.994915,
                                                 21.981943,                         -36.021129,                         31.989829};
    for (int k = 0; k < 200; k++)
    {
      const double shift = -300.0 + k * 3.0157;
      if (std::any_of(neighbours.begin(), neighbours.end(), [&](double n) { return std::abs(shift - n) < 0.3; })) continue;
      if (count_pairs(precursor_mass + shift) >= observed) return false;
    }
    return true;
  }

  void FLASHExtenderAlgorithm::calculatePrecursorMass_(const ProteinHit& hit, const std::map<int, std::vector<Size>>& best_path_map, HitInformation& hi)
  {
    hi.calculated_precursor_mass_ = -1;
    if (hi.protein_start_position_ < 0 || hi.protein_end_position_ < 0) return;
    // the inference needs a path from each terminus. The check above already implies both (the start position is set only by a prefix
    // path, the end position only by a suffix path); this keeps at() from throwing inside the parallel loop should that ever change
    if (best_path_map.find(0) == best_path_map.end() || best_path_map.find(1) == best_path_map.end()) return;

    const auto& bp0 = best_path_map.at(0);
    const auto& bp1 = best_path_map.at(1);

    // PeakGroup pg;
    // pg.setQscore(.8);
    const int terminal_score_threshold = 1;   //(int)round(FLASHTaggerAlgorithm::getNodeScore(pg) * 3); // n term or c term should have at least this
                                              //score to be considered = three high quality masses
    int max_score = terminal_score_threshold; // need to exceed 20 for precursor correction
    double min_tol = 1000;
    Size pro_size = hi.pro_mass_map_[0].size();
    int min_excessive_aa = pro_size;
    for (auto iter0 = bp0.rbegin(); iter0 != bp0.rend(); iter0++)
    {
      int proi = getProIndex_(*iter0, pro_size);
      int proj = 0;
      auto iter1 = bp1.rbegin();
      for (; iter1 != bp1.rend(); iter1++)
      {
        proj = getProIndex_(*iter1, pro_size);
        if (proi + proj >= (int)hit.getSequence().size()) break;
      }
      if (iter1 == bp1.rend()) continue; // no bp1 vertex complements this bp0 vertex to full coverage
      if (getScore_(*iter0) < terminal_score_threshold || getScore_(*iter1) < terminal_score_threshold) continue; //
      int excessive_aa = (proi + proj) - (int)hit.getSequence().size();

      if (excessive_aa < 0) continue;
      if (excessive_aa > min_excessive_aa) continue;
      if (excessive_aa < min_excessive_aa)
      {
        max_score = terminal_score_threshold;
        min_tol = 1000;
        min_excessive_aa = excessive_aa;
      }

      int score = getScore_(*iter0) + getScore_(*iter1);
      if (max_score > score) continue;

      if (max_score == score)
      {
        double tol = std::max(hi.tol_spec_map_.at(0)[getNodeIndex_(*iter0, pro_size)].getIntensity(),
                              hi.tol_spec_map_.at(1)[getNodeIndex_(*iter1, pro_size)].getIntensity());
        if (tol > min_tol) continue;
        min_tol = tol;
      }
      max_score = score;
      hi.calculated_precursor_mass_
        = hi.node_spec_map_.at(0)[getNodeIndex_(*iter0, pro_size)].getMZ() + hi.node_spec_map_.at(1)[getNodeIndex_(*iter1, pro_size)].getMZ();

      if (proi + proj > (int)hit.getSequence().size())
      {
        hi.calculated_precursor_mass_ += hi.pro_mass_map_[0].back() - (hi.pro_mass_map_[0][proi] + hi.pro_mass_map_[1][proj]);
      }
    }
    if (hi.calculated_precursor_mass_ > 0) hi.calculated_precursor_mass_ += i2f_mass;
  }

  void FLASHExtenderAlgorithm::getProMasses_(const ProteinHit& hit, std::vector<double>& pro_masses, int mode)
  {
    pro_masses.reserve(hit.getSequence().size() + 1);
    pro_masses.push_back(0);

    auto seq = hit.getSequence();

    if (mode == 0) seq = seq.reverse();
    for (const auto& aa : seq)
    {
      if (aa == 'X') pro_masses.push_back(pro_masses.back()); // repeat the previous mass
      else
        pro_masses.push_back(pro_masses.back() + ResidueDB::getInstance()->getResidue(aa)->getMonoWeight(Residue::Internal));
    }
  }

  double FLASHExtenderAlgorithm::getSpecMassSpan_(const std::vector<Size>& path, const MSSpectrum& node_spec, int pro_mass_size) const
  {
    return node_spec[(getNodeIndex_(path[0], pro_mass_size))].getMZ();
  }

  double FLASHExtenderAlgorithm::getProteinMassSpan_(const std::vector<Size>& path, const std::vector<double>& pro_masses) const
  {
    double minmass = 0, maxmass = 0;
    maxmass = pro_masses[getProIndex_(path[0], pro_masses.size())];

    for (int i = path.size() - 1; i >= 0; i--)
    {
      if (getNodeIndex_(path[i], pro_masses.size()) != 0) break;
      minmass = pro_masses[getProIndex_(path[i], pro_masses.size())];
    }
    return std::abs(maxmass - minmass);
  }

  int FLASHExtenderAlgorithm::getModifiedAACount_(const std::vector<Size>& path) const
  {
    int cntr = 0;
    for (auto vertex : path)
    {
      if (getBlindModNumber_(vertex) == 0) continue;
      cntr++;
    }
    return cntr;
  }

  int FLASHExtenderAlgorithm::getProteinLength_(const std::vector<Size>& path, const std::vector<double>& pro_masses) const
  {
    int minindex = 0, maxindex = 0;
    maxindex = getProIndex_(path[0], pro_masses.size());

    for (int i = path.size() - 1; i >= 0; i--)
    {
      if (getNodeIndex_(path[i], pro_masses.size()) != 0) break;
      minindex = getProIndex_(path[i], pro_masses.size());
    }
    return std::abs(maxindex - minindex);
  }


  int FLASHExtenderAlgorithm::getBestMassShift_(const std::vector<int>& spec_vec,
                                                const std::unordered_set<int>& pro_vec)
  {
    int ret = 0;
    int match_cntr = 0;
    if (spec_vec.empty() || pro_vec.empty()) return ret;

    // counts per integer shift in a reused array instead of a hash map; the visiting order (and so the winner among equal counts) is
    // unchanged
    const auto [smin, smax] = std::minmax_element(spec_vec.begin(), spec_vec.end());
    int pmin = *pro_vec.begin(), pmax = pmin;
    for (int p : pro_vec)
    {
      pmin = std::min(pmin, p);
      pmax = std::max(pmax, p);
    }
    const int offset = pmin - *smax;
    const Size range = (Size)(pmax - *smin - offset + 1);
    static thread_local std::vector<int> deltas;
    if (deltas.size() < range) deltas.assign(range, 0);
    for (int s : spec_vec)
    {
      for (int p : pro_vec)
      {
        const int d = p - s;
        const int cnt = ++deltas[d - offset];
        if (cnt > match_cntr)
        {
          match_cntr = cnt;
          ret = d;
        }
      }
    }
    for (int s : spec_vec) // reset only the touched entries
    {
      for (int p : pro_vec) deltas[p - s - offset] = 0;
    }
    return ret;
  }


  void FLASHExtenderAlgorithm::updateHitInformation_(const DeconvolvedSpectrum& dspec, HitInformation& hi, double max_mass)
  {
    // 0 for suffix 1 for prefix 2 for suffix and prefix if precursor mass is available
    double max_proteoform_mass = -1;
    float max_score = -(float)max_path_score_;

    std::vector<int> isos {0};
    if (hi.mode_ == 2 && hi.calculated_precursor_mass_ > 0 && ! hi.precursor_as_given_)
    {
      for (int i = 1; i <= allowed_isotope_error_; ++i)
      {
        isos.push_back(i);
        isos.push_back(-i);
      }
    }
    for (int iso : isos)
    {
      double proteoform_mass = hi.calculated_precursor_mass_ + iso * Constants::C13C12_MASSDIFF_U;
      // (mass, score, tolerance, series) of every fragment mass of this mode, sorted together so that the series stays with its mass
      struct Node
      {
        double mass;
        float score;
        double tol;
        int series;
      };
      std::vector<Node> nodes;
      nodes.reserve(dspec.size() * ion_types_str_.size() + 1);

      for (const auto& pg : dspec)
      {
        const auto& shifts = (hi.mode_ == 0)   ? suffix_shifts_
                             : (hi.mode_ == 1) ? prefix_shifts_
                             : (hi.mode_ == 2) ? prefix_shifts_
                                               : std::vector<double> {};
        const auto& series = hi.mode_ == 0 ? suffix_series_ : prefix_series_;

        for (Size si = 0; si < shifts.size(); si++)
        {
          double mass = pg.getMonoMass() - shifts[si];
          if (mass <= 0 || mass > max_mass + 1) continue;
          nodes.push_back({mass, (float)std::max(1, FLASHTaggerAlgorithm::getNodeScore(pg)), tol_ * pg.getMonoMass(), series[si]});
        }

        if (hi.mode_ == 2 && hi.calculated_precursor_mass_ > 0)
        {
          for (Size si = 0; si < suffix_shifts_.size(); si++)
          {
            double mass = pg.getMonoMass() - suffix_shifts_[si];
            if (mass <= 0 || mass >= proteoform_mass - i2f_mass) continue;
            double prefix_mass = proteoform_mass - i2f_mass - mass;
            // the converted mass carries the error of the precursor mass as well as that of the fragment mass (both at the input tolerance)
            nodes.push_back({prefix_mass, (float)std::max(1, FLASHTaggerAlgorithm::getNodeScore(pg)), tol_ * (pg.getMonoMass() + proteoform_mass),
                             suffix_series_[si]});
          }
        }
      }
      std::stable_sort(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) { return a.mass < b.mass; });
      MSSpectrum node_spec, tol_spec;
      node_spec.reserve(nodes.size());
      tol_spec.reserve(nodes.size());
      for (const auto& n : nodes)
      {
        node_spec.emplace_back(n.mass, n.score);
        tol_spec.emplace_back(n.mass, n.tol);
      }

      MSSpectrum t_node_spec, t_tol_spec;
      std::vector<int> t_series;
      t_node_spec.reserve(node_spec.size() + 1);
      t_tol_spec.reserve(node_spec.size() + 1);
      t_series.reserve(node_spec.size() + 1);

      t_node_spec.emplace_back(0, 0);
      t_tol_spec.emplace_back(0, 0);
      t_series.push_back(-1);
      float overlapped_score = 0;

      for (Size k = 0; k < node_spec.size(); k++)
      {
        const auto& p = node_spec[k];
        double mass = p.getMZ();
        float score = p.getIntensity();
        int series = nodes[k].series;
        // the half-width of the node merged so far (t_node_spec.back()), not that of the raw mass k - 1: after a merge they differ, and
        // comparing or merging with the raw one missed overlaps and shifted the merged centre. At k == 0, back() is the (0, 0) entry.
        const double prev_margin = t_tol_spec.back().getIntensity();
        double margin = tol_spec[k].getIntensity();

        if (mass - margin < t_node_spec.back().getMZ() + prev_margin) // they are the same
        {
          float prev_score = t_node_spec.back().getIntensity();

          const double orig_margin = margin;
          margin = (mass + orig_margin - t_node_spec.back().getMZ() + prev_margin) / 2.0;
          mass = (mass + orig_margin + t_node_spec.back().getMZ() - prev_margin) / 2.0;
          if (t_node_spec.size() > 1)
          {
            if (t_series.back() != series) series = -1; // a merged b/y and c/z mass: no series
            t_node_spec.pop_back();
            t_tol_spec.pop_back();
            t_series.pop_back();
          }
          score = multi_ion_score + std::max(prev_score, score);
          // score = suffix_prefix == prev_suffix_prefix ? multi_ion_score + std::max(prev_score, score) : multi_ion_score + prev_score + score;
          overlapped_score += score;
        }

        if (mass <= 0) continue;

        t_node_spec.emplace_back(mass, score);
        t_tol_spec.emplace_back(mass, margin);
        t_series.push_back(series);
      }

      if (max_score < overlapped_score - multi_ion_score)
      {
        hi.node_spec_map_[hi.mode_] = t_node_spec;
        hi.tol_spec_map_[hi.mode_] = t_tol_spec;
        hi.node_series_map_[hi.mode_] = t_series;
        max_score = overlapped_score;
        max_proteoform_mass = proteoform_mass;
      }
    }

    hi.node_spec_map_[hi.mode_].sortByPosition();
    hi.tol_spec_map_[hi.mode_].sortByPosition();
    if (hi.calculated_precursor_mass_ > 0)
    {
      hi.calculated_precursor_mass_ = max_proteoform_mass;
      hi.node_spec_map_[hi.mode_].emplace_back(hi.calculated_precursor_mass_ - i2f_mass, 1);
      hi.tol_spec_map_[hi.mode_].emplace_back(hi.calculated_precursor_mass_ - i2f_mass, tol_ * hi.calculated_precursor_mass_);
      hi.node_series_map_[hi.mode_].push_back(-1);
    }
  }

  void FLASHExtenderAlgorithm::run_(const ProteinHit& hit,
                                    HitInformation& hi,
                                    const std::vector<FLASHHelperClasses::Tag>& matched_tags,
                                    std::map<int, std::vector<Size>>& all_paths_per_mode,
                                    int max_mod_cntr_for_last_mode) // per hit
  {
    std::vector<std::vector<int>> tag_edges(4); // pro start end indices and node start end indices
    std::set<Size> sinks;
    std::set<String> upper_case_seqs;

    for (const auto& tag : matched_tags)
    {
      if (std::any_of(tag.getSequence().begin(), tag.getSequence().end(), ::islower)) continue;
      upper_case_seqs.insert(tag.getSequence());
    }

    if (! matched_tags.empty())
    {
      for (auto& edge : tag_edges)
      {
        edge.reserve(matched_tags.size() * 2);
      }
    }
    const auto& node_spec = hi.node_spec_map_[hi.mode_];
    const auto& pro_masses = hi.pro_mass_map_[hi.mode_];
    // cache current-mode map entries so the recursion does not re-resolve these std::map lookups on every frame (entries are stable for this mode)
    hi.cur_node_spec_ = &node_spec;
    hi.cur_pro_masses_ = &pro_masses;
    hi.cur_tol_spec_ = &hi.tol_spec_map_[hi.mode_];
    hi.cur_oriented_seq_ = hit.getSequence();
    if (hi.mode_ == 0) hi.cur_oriented_seq_.reverse();
    hi.cur_node_series_ = &hi.node_series_map_[hi.mode_];
    {
      // weight of a matched mass at protein index k (the bond after k residues in this mode's orientation)
      const String& os = hi.cur_oriented_seq_;
      const int L = (int)os.size();
      for (int sr = 0; sr < 2; sr++)
      {
        auto& w = hi.cur_bond_weight_[sr];
        w.assign(L + 1, 0);
        if (cleavage_weight_ <= 0) continue;
        for (int k = 1; k < L; k++)
        {
          const char before = hi.mode_ == 0 ? os[k] : os[k - 1]; // N-terminal side of the bond
          const char after = hi.mode_ == 0 ? os[k - 1] : os[k];
          if (before < 'A' || before > 'Z' || after < 'A' || after > 'Z') continue;
          const double v = cleavage_weight_ * (cleavage_before_[sr][before - 'A'] + cleavage_after_[sr][after - 'A']);
          w[k] = std::max(-4, std::min(4, (int)std::lround(v)));
        }
      }
    }
    hi.dag_ = FLASHHelperClasses::DAG((1 + max_var_mod_cntr_) * (1 + node_spec.size()) * (1 + pro_masses.size()) * (1 + max_blind_mod_cntr_)
                                      * (1 + max_path_score_ - min_path_score_));

    bool tag_found = false;
    const String& seq = hit.getSequence();

    for (const auto& tag : matched_tags)
    {
      if (hi.mode_ == 2 && hi.calculated_precursor_mass_ <= 0) continue;
      const auto& upper_tag_seq = tag.getUppercaseSequence();
      bool has_lower_case = std::any_of(tag.getSequence().begin(), tag.getSequence().end(), ::islower);
      if (has_lower_case && upper_case_seqs.find(upper_tag_seq) != upper_case_seqs.end()) continue;

      tag_found = true;
      std::vector<int> positions;
      std::vector<double> masses;

      FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(positions, masses, -1, seq, tag);
      auto tag_masses = tag.getMzs();
      std::sort(tag_masses.begin(), tag_masses.end());

      double seq_mass = 0;
      if (has_lower_case) { seq_mass = AASequence::fromString(upper_tag_seq, true).getMonoWeight(Residue::Internal); }
      std::vector<double> start_masses, end_masses, start_tols, end_tols;
      if (tag.getCtermMass() >= 0) // suffix
      {
        for (const auto& shift : suffix_shifts_)
        {
          double start_mass = tag_masses[0] - shift;
          double end_mass = tag_masses.back() - shift;

          if (seq_mass > 0)
          {
            if (std::abs(start_mass + seq_mass - end_mass) > end_mass * tol_) end_mass = start_mass + seq_mass;
          }

          // converted with the precursor mass (mode 2), a suffix mass also carries the precursor mass error
          const double converted = hi.mode_ == 2 ? hi.calculated_precursor_mass_ : 0;
          start_tols.push_back((start_mass + converted) * tol_);
          end_tols.push_back((end_mass + converted) * tol_);

          if (hi.mode_ == 2)
          {
            start_mass = hi.calculated_precursor_mass_ - i2f_mass - start_mass;
            end_mass = hi.calculated_precursor_mass_ - i2f_mass - end_mass;
            start_masses.push_back(end_mass);
            end_masses.push_back(start_mass);
          }
          else
          {
            start_masses.push_back(start_mass);
            end_masses.push_back(end_mass);
          }
        }
      }
      else // prefix
      {
        for (const auto& shift : prefix_shifts_)
        {
          double start_mass = tag_masses[0] - shift;
          double end_mass = tag_masses.back() - shift;

          if (seq_mass > 0)
          {
            if (std::abs(start_mass + seq_mass - end_mass) > end_mass * tol_) end_mass = start_mass + seq_mass;
          }

          start_tols.push_back(start_mass * tol_);
          end_tols.push_back(end_mass * tol_);
          start_masses.push_back(start_mass);
          end_masses.push_back(end_mass);
        }
      }

      bool same_terminal_tag = hi.mode_ == 2 || (hi.mode_ == 0 && tag.getCtermMass() >= 0) || (hi.mode_ == 1 && tag.getNtermMass() >= 0);

      for (Size l = 0; l < start_masses.size(); l++)
      {
        int highest_score_start = -1, highest_score_end = -1;
        if (same_terminal_tag)
        {
          double delta_start = start_tols[l];
          double delta_end = end_tols[l];

          highest_score_start = node_spec.findHighestInWindow(start_masses[l], delta_start, delta_start);
          highest_score_end = node_spec.findHighestInWindow(end_masses[l], delta_end, delta_end);

          if (highest_score_start < 0 || highest_score_end < 0 || highest_score_start >= (int)node_spec.size()
              || highest_score_end >= (int)node_spec.size())
            continue;
        }
        for (int pos : positions)
        {
          if (hi.mode_ == 0) // suffix inverted
          {
            pos = (int)pro_masses.size() - 1 - pos; // invert pos
            if (pos - tag.getLength() >= 0 && pos < (int)pro_masses.size())
            {
              tag_edges[0].emplace_back(pos - tag.getLength());
              tag_edges[1].emplace_back(pos);
              tag_edges[2].emplace_back(highest_score_start);
              tag_edges[3].emplace_back(highest_score_end); // check...
            }
          }
          else
          {
            if (pos >= 0 && pos + tag.getLength() < pro_masses.size())
            {
              tag_edges[0].emplace_back(pos);
              tag_edges[1].emplace_back(pos + tag.getLength()); // this can be much faster...
              tag_edges[2].emplace_back(highest_score_start);
              tag_edges[3].emplace_back(highest_score_end);
            }
          }
        }
      }
    }
    // Tags of the other terminus cannot be placed on this mode's spectrum (node index -1). Such an edge only bounds a free extension
    // from the terminus, and findPathsAlongTagEndPoints replaces all sinks with that extension, so mixing them with this mode's own
    // tags discarded the paths through those tags. Use this mode's own tags first; only if they give no path with a positive score
    // (e.g. all of them sit at wrong positions), fall back to the free extension bounded by the other terminus' tags.
    std::vector<std::vector<int>> own_tag_edges(tag_edges.size()), other_tag_edges(tag_edges.size());
    for (Size e = 0; e < tag_edges[2].size(); e++)
    {
      auto& edges = tag_edges[2][e] >= 0 ? own_tag_edges : other_tag_edges;
      for (Size k = 0; k < tag_edges.size(); k++) edges[k].push_back(tag_edges[k][e]);
    }
    if (own_tag_edges[2].empty() || other_tag_edges[2].empty()) { constructDAG_(sinks, hi, tag_edges, max_mod_cntr_for_last_mode, tag_found); }
    else
    {
      constructDAG_(sinks, hi, own_tag_edges, max_mod_cntr_for_last_mode, tag_found);
      if (std::none_of(sinks.begin(), sinks.end(), [this](Size v) { return getScore_(v) > 0; }))
      {
        sinks.clear();
        hi.dag_ = FLASHHelperClasses::DAG(hi.dag_.size());
        hi.visited_.clear();
        constructDAG_(sinks, hi, other_tag_edges, max_mod_cntr_for_last_mode, tag_found);
      }
    }
    Size src = getVertex_(0, 0, 0, 0, 0, pro_masses.size());
    std::vector<int> max_scores(max_blind_mod_cntr_ + 1, 0);
    for (Size sink : sinks)
    {
      int num_mod = getBlindModNumber_(sink);
      if (sink == src || getScore_(sink) < max_scores[num_mod]) continue; // getBlindModNumber_(sink) == 0 ||
      if (hi.calculated_precursor_mass_ > 0 && getNodeIndex_(sink, pro_masses.size()) < (int)node_spec.size() - 1) continue;
      max_scores[num_mod] = getScore_(sink);
    }
    std::vector<std::vector<std::vector<Size>>> paths(max_blind_mod_cntr_ + 1, std::vector<std::vector<Size>>());

    for (Size sink : sinks)
    {
      int num_mod = getBlindModNumber_(sink);
      if (sink == src || getScore_(sink) < max_scores[num_mod]) continue; //
      if (hi.calculated_precursor_mass_ > 0 && getNodeIndex_(sink, pro_masses.size()) < (int)node_spec.size() - 1) continue;
      // std::vector<std::vector<Size>> sub_paths;
      hi.dag_.findAllPaths(sink, src, paths[num_mod], 0);
    }
    for (int num_mod = 0; num_mod <= max_blind_mod_cntr_; num_mod++)
    {
      for (const auto& path : paths[num_mod])
      {
        double mass = getSpecMassSpan_(path, node_spec, pro_masses.size());
        double pro_mass = getProteinMassSpan_(path, pro_masses);
        int pro_len = getProteinLength_(path, pro_masses);

        auto iter = all_paths_per_mode.find(num_mod);

        if (iter == all_paths_per_mode.end()) { all_paths_per_mode[num_mod] = path; }
        else
        {
          double mod_mass_
            = num_mod == 0 ? 0 : std::abs(getSpecMassSpan_(iter->second, node_spec, pro_masses.size()) - getProteinMassSpan_(iter->second, pro_masses));
          double mod_mass = num_mod == 0 ? 0 : std::abs(mass - pro_mass);
          if (mod_mass_ < mod_mass) continue;

          // Prefer the path with smaller mod mass
          if (mod_mass_ > mod_mass
              || (mod_mass_ == mod_mass
                  && (getProteinLength_(iter->second, pro_masses) < pro_len
                      || (getProteinLength_(iter->second, pro_masses) == pro_len && getModifiedAACount_(path) < getModifiedAACount_(iter->second)))))
          {
            all_paths_per_mode[num_mod] = path;
          }
        }
      }
    }
  }

  void FLASHExtenderAlgorithm::run(std::vector<ProteinHit>& hits,
                                   const DeconvolvedSpectrum& dspec, const std::vector<int>& spec_vec, const std::vector<std::unordered_set<int>>& vec_pro,
                                   const std::vector<std::unordered_set<int>>& rev_vec_pro, const std::vector<FLASHHelperClasses::Tag>& tags, double ppm, bool multiple_hits_per_spec)
  {
    if (hits.empty()) return;
    // setLogType(CMD);
    //
    //var_mods_[27.994915	] = ResidueModification();
    ion_types_str_ = param_.getValue("ion_type").toStringVector();
    std::sort(ion_types_str_.begin(), ion_types_str_.end());
    for (const auto& ion_str : ion_types_str_)
    {
      if (ion_str == "a") { prefix_shifts_.push_back(Residue::getInternalToAIon().getMonoWeight()); prefix_series_.push_back(0); }
      else if (ion_str == "b") { prefix_shifts_.push_back(Residue::getInternalToBIon().getMonoWeight()); prefix_series_.push_back(0); }
      else if (ion_str == "c") { prefix_shifts_.push_back(Residue::getInternalToCIon().getMonoWeight()); prefix_series_.push_back(1); }
      else if (ion_str == "x") { suffix_shifts_.push_back(Residue::getInternalToXIon().getMonoWeight()); suffix_series_.push_back(0); }
      else if (ion_str == "y") { suffix_shifts_.push_back(Residue::getInternalToYIon().getMonoWeight()); suffix_series_.push_back(0); }
      else if (ion_str == "z") { suffix_shifts_.push_back(Residue::getInternalToZIon().getMonoWeight()); suffix_series_.push_back(1); }
      else if (ion_str == "zp1") { suffix_shifts_.push_back(Residue::getInternalToZp1Ion().getMonoWeight()); suffix_series_.push_back(1); }
      else if (ion_str == "zp2") { suffix_shifts_.push_back(Residue::getInternalToZp2Ion().getMonoWeight()); suffix_series_.push_back(1); }
      else { continue; }
    }

    tol_ = ppm / 1e6;
    proteoform_hits_.clear();

    std::vector<double> mzs;
    std::vector<int> scores;

    proteoform_hits_.reserve(hits.size());
    if (dspec.getPrecursorPeakGroup().getMonoMass() > 0) { given_precursor_mass_ = dspec.getPrecursorPeakGroup().getMonoMass(); }

    // meta values set or read by name look the name up in the global registry under a critical section; use indices in the loop
    auto& registry = MetaInfoInterface::metaRegistry();
    const UInt mi_tag_indices = registry.registerName("TagIndices");
    const UInt mi_fasta_index = registry.registerName("FastaIndex");
    const UInt mi_modification_ids = registry.registerName("ModificationIDs");
    const UInt mi_modification_accs = registry.registerName("ModificationACCs");
    const UInt mi_modifications = registry.registerName("Modifications");
    const UInt mi_modification_starts = registry.registerName("ModificationStarts");
    const UInt mi_modification_ends = registry.registerName("ModificationEnds");
    const UInt mi_blind_mod_count = registry.registerName("BlindModCount");
    const UInt mi_variable_modifications = registry.registerName("VariableModifications");
    const UInt mi_modification_variable = registry.registerName("ModificationVariable");
    const UInt mi_matched_aa = registry.registerName("MatchedAA");
    const UInt mi_start_position = registry.registerName("StartPosition");
    const UInt mi_end_position = registry.registerName("EndPosition");
    const UInt mi_given_mass = registry.registerName("GivenMass");
    const UInt mi_mass = registry.registerName("Mass");
    const UInt mi_rt = registry.registerName("RT");
    const UInt mi_num_mass = registry.registerName("NumMass");
    const UInt mi_precursor_score = registry.registerName("PrecursorScore");
    const UInt mi_precursor_snr = registry.registerName("PrecursorSNR");
    const UInt mi_mass_by_fragments = registry.registerName("ProteoformMassByFragmentMass");

    startProgress(0, (int)hits.size(), "running FLASHExtender ...");
    std::vector<std::pair<int, ProteinHit>> found; // (candidate index, proteoform hit)

    std::vector<double> sorted_masses; // for the complementary-pair check of inferred precursor masses
    sorted_masses.reserve(dspec.size());
    for (const auto& pg : dspec) sorted_masses.push_back(pg.getMonoMass());
    std::sort(sorted_masses.begin(), sorted_masses.end());

  #pragma omp parallel for default(none) schedule(dynamic) if (! omp_in_parallel()) shared(hits, tags, dspec, multiple_hits_per_spec, i2f_mass, spec_vec, sorted_masses, vec_pro, rev_vec_pro, std::cout, found, mi_tag_indices, mi_fasta_index, mi_modification_ids, mi_modification_accs, mi_modifications, mi_modification_starts, mi_modification_ends, mi_blind_mod_count, mi_variable_modifications, mi_modification_variable, mi_matched_aa, mi_start_position, mi_end_position, mi_given_mass, mi_mass, mi_rt, mi_num_mass, mi_precursor_score, mi_precursor_snr, mi_mass_by_fragments)
    for (int i = 0; i < ((int)hits.size()); i++)
    {
      nextProgress();
      auto& hit = hits[i];
      HitInformation hi;
      int total_score = 0;
      std::vector<int> mod_starts, mod_ends;
      std::vector<double> mod_masses, mod_tols;
      std::vector<int> mod_vars;                 // 1: the modification is a variable one (named in mod_var_names), 0: blind
      std::vector<String> mod_var_names, mod_var_accs;
      int max_nterm_index = 0, max_cterm_rindex = 0;

      std::map<int, std::map<int, std::vector<Size>>> all_path_map; // mode, num_mod, path
      std::map<int, std::vector<Size>> best_path_map;               // mode, best paths

      std::vector<int> used_mode;

      if (! hit.metaValueExists(mi_tag_indices)) continue;

      const std::vector<int>& tag_indices = hit.getMetaValue(mi_tag_indices).toIntList();

      std::vector<FLASHHelperClasses::Tag> matched_tags;
      for (const auto ti : tag_indices)
      {
        matched_tags.push_back(tags[ti]);
      }

      std::map<int, std::set<int>> matched_position_map;
      bool precursor_by_fragment = false;

      hi.n_best_shift_ = getBestMassShift_(spec_vec, vec_pro[hit.getMetaValue(mi_fasta_index)]);
      hi.c_best_shift_ = getBestMassShift_(spec_vec, rev_vec_pro[hit.getMetaValue(mi_fasta_index)]);

      // The result of the terminal modes (0: suffix, 1: prefix), kept when mode 2 starts. Mode 2 (both termini, constrained by the
      // precursor mass) replaces their score; if it finds no path, or only one that scores <= 0 (e.g. a fragment truncated at both
      // ends that the precursor mass does not admit, or a wrong precursor mass), the terminal result is restored with an
      // undetermined mass instead of dropping the hit.
      bool terminal_kept = false;
      int terminal_score = 0;
      std::vector<int> terminal_used_mode, terminal_mod_starts, terminal_mod_ends, terminal_mod_vars;
      std::vector<double> terminal_mod_masses, terminal_mod_tols;
      std::vector<String> terminal_mod_var_names, terminal_mod_var_accs;
      std::map<int, std::set<int>> terminal_matched_position_map;
      int terminal_start = -1, terminal_end = -1;
      // suffix fragment masses matched on the suffix path before its first shift, i.e. measured from where that path starts
      int c_ladder = 0, terminal_c_ladder = 0;
      bool retried_with_given_mass = false;

      for (hi.mode_ = 0; hi.mode_ <= 2; hi.mode_++)
      {
        int max_mod_cntr_for_last_mode = -1;
        if (hi.mode_ == 2)
        {
          terminal_kept = true;
          terminal_score = total_score;
          terminal_used_mode = used_mode;
          terminal_mod_starts = mod_starts;
          terminal_mod_ends = mod_ends;
          terminal_mod_masses = mod_masses;
          terminal_mod_tols = mod_tols;
          terminal_mod_vars = mod_vars;
          terminal_mod_var_names = mod_var_names;
          terminal_mod_var_accs = mod_var_accs;
          terminal_matched_position_map = matched_position_map;
          terminal_start = hi.protein_start_position_;
          terminal_end = hi.protein_end_position_;
          terminal_c_ladder = c_ladder;
        }
        if (hi.mode_ == 2 && hi.calculated_precursor_mass_ <= 0)
        { // const ProteinHit& hit,
          if (! skip_precursor_inference_ && ! retried_with_given_mass && max_nterm_index + max_cterm_rindex >= (int)hit.getSequence().size())
            calculatePrecursorMass_(hit, best_path_map, hi);
          max_mod_cntr_for_last_mode = std::min(max_blind_mod_cntr_, (int)std::count(mod_vars.begin(), mod_vars.end(), 0) + 1);

          if (hi.calculated_precursor_mass_ <= 0) hi.calculated_precursor_mass_ = given_precursor_mass_;
          else
            precursor_by_fragment = true;
          if (hi.calculated_precursor_mass_ <= 0) break;
        }

        auto& pro_masses = hi.pro_mass_map_[hi.mode_] = std::vector<double>();
        // auto& node_spec = hi.node_spec_map_[hi.mode_] = MSSpectrum();
        // auto& tol_spec = hi.tol_spec_map_[hi.mode_] = MSSpectrum();

        getProMasses_(hit, pro_masses, hi.mode_);
        updateHitInformation_(dspec, hi, pro_masses.back());

        hi.visited_.clear(); // dag_ is rebuilt per mode

        run_(hit, hi, matched_tags, all_path_map[hi.mode_], max_mod_cntr_for_last_mode);

        if (hi.mode_ < 2)
        {
          const auto paths_c = all_path_map.find(0);
          const auto paths_n = all_path_map.find(1);

          std::map<int, int> nscores, cscores; // mod and score

          if (paths_n != all_path_map.end())
          {
            for (const auto& [mod, path] : paths_n->second)
            {
              if (path.empty()) continue;
              nscores[mod] = getScore_(path[0]);
            }
          }
          if (paths_c != all_path_map.end())
          {
            for (const auto& [mod, path] : paths_c->second)
            {
              if (path.empty()) continue;
              cscores[mod] = getScore_(path[0]);
            }
          }
          int max_score = 0;

          if (nscores.empty()) // only c term
          {
            for (const auto& [mod, path] : paths_c->second)
            {
              if (path.empty()) continue;
              if (max_score >= getScore_(path[0])) continue;
              max_score = getScore_(path[0]);
              best_path_map[0] = path;
            }
          }
          else if (cscores.empty()) // only n term
          {
            for (const auto& [mod, path] : paths_n->second)
            {
              if (path.empty()) continue;
              if (max_score >= getScore_(path[0])) continue;
              max_score = getScore_(path[0]);
              best_path_map[1] = path;
            }
          }
          else // both terms
          {
            for (int mc = 0; mc <= max_blind_mod_cntr_; mc++)
            {
              if (cscores.find(mc) == cscores.end()) continue;
              const auto& cpath = paths_c->second[mc];
              for (int mn = 0; mc + mn <= max_blind_mod_cntr_; mn++)
              {
                if (nscores.find(mn) == nscores.end()) continue;
                int sum_score = nscores[mn] + cscores[mc];
                if (max_score >= sum_score) continue;
                max_score = sum_score;
                best_path_map[0] = cpath;
                best_path_map[1] = paths_n->second[mn];
              }
            }
          }
        }
        else if (hi.mode_ == 2)
        {
          int max_score = 0;
          const auto paths = all_path_map.find(2);

          for (const auto& [mod, path] : paths->second)
          {
            if (path.empty()) continue;
            if (max_score >= getScore_(path[0])) continue;
            max_score = getScore_(path[0]);
            best_path_map[2] = path;
          }
          mod_starts.clear();
          mod_ends.clear();
          mod_masses.clear();
          mod_tols.clear();
          mod_vars.clear();
          mod_var_names.clear();
          mod_var_accs.clear();
          //
          //        int mode0_score = best_path_map.find(0) == best_path_map.end() ||  best_path_map[0].empty() ? 0 :  getScore_(best_path_map[0][0]);
          //        int mode1_score = best_path_map.find(1) == best_path_map.end() ||  best_path_map[1].empty() ? 0 :  getScore_(best_path_map[1][0]);
          //
          //        if (max_score < mode0_score + mode1_score)
          //        {
          //          precursor_by_fragment = false;
          //          hi.calculated_precursor_mass_ = -1;
          //          break;
          //        }
        }
        if (hi.mode_ == 0) continue;
        total_score = 0;
        // find the best paths per mode. Mode 0 and 1 should be considered together (since the modification counts for N C term paths should be summed
        // up).
        // number of modifications recorded before each mode's path was read, so that the modifications of a mode dropped below can be
        // removed again
        std::map<int, Size> mods_before;
        auto erase_mods = [&](Size from, Size to) {
          mod_starts.erase(mod_starts.begin() + from, mod_starts.begin() + to);
          mod_ends.erase(mod_ends.begin() + from, mod_ends.begin() + to);
          mod_masses.erase(mod_masses.begin() + from, mod_masses.begin() + to);
          mod_tols.erase(mod_tols.begin() + from, mod_tols.begin() + to);
          mod_vars.erase(mod_vars.begin() + from, mod_vars.begin() + to);
          mod_var_names.erase(mod_var_names.begin() + from, mod_var_names.begin() + to);
          mod_var_accs.erase(mod_var_accs.begin() + from, mod_var_accs.begin() + to);
        };
        for (int m = hi.mode_ == 2 ? 2 : 0; m <= hi.mode_; m++)
        {
          if (best_path_map.empty() || best_path_map.find(m) == best_path_map.end() || best_path_map[m].empty()) continue;
          mods_before[m] = mod_starts.size();
          if (m == 0) c_ladder = 0;
          auto& best_path = best_path_map[m];
          auto& t_pro_masses = hi.pro_mass_map_[m];
          auto& t_node_spec = hi.node_spec_map_[m];

          double prev_mass_shift = 0;
          int prev_mod_count = 0, prev_var_count = 0;
          String oriented_seq = hit.getSequence();
          if (m == 0) oriented_seq.reverse();
          int pre_pro_index = 0;
          int pre_node_index = 0;
          double mod_mass = 0; //, total_mod_mass = 0;
          // int total_mod_count = getBlindModNumber_(*best_path.begin());

          for (auto iter = best_path.rbegin(); iter != best_path.rend(); iter++)
          {
            auto pro_index = getProIndex_(*iter, t_pro_masses.size());
            auto node_index = getNodeIndex_(*iter, t_pro_masses.size());
            auto mass_shift = t_node_spec[node_index].getMZ() - t_pro_masses[pro_index];
            const int var_count = getVarModNumber_(*iter);
            auto mod_count = getBlindModNumber_(*iter) + var_count; // blind and variable modifications so far

            if (node_index == 0)
            {
              if (m > 0) hi.protein_start_position_ = pro_index;
              if (m == 0) hi.protein_end_position_ = (int)hit.getSequence().size() - pro_index; //
              if (mod_count == 0) prev_mass_shift = mass_shift;
            }

            int pro_seq_index = m > 0 ? pro_index : ((int)hit.getSequence().size() - pro_index);
            if (node_index > 0 && t_node_spec[node_index].getMZ() != hi.calculated_precursor_mass_ - i2f_mass)
            {
              matched_position_map[m].insert(pro_seq_index);
              if (m == 0 && mod_count == 0) c_ladder++;
            }

            if (m == 0) max_cterm_rindex = std::max(max_cterm_rindex, pro_index);
            if (m == 1) max_nterm_index = std::max(max_nterm_index, pro_index);
            if (m == 2) hi.protein_end_position_ = pro_index;
            if (mod_count != prev_mod_count)
            {
              mod_mass = mass_shift - prev_mass_shift;
              // total_mod_mass += mod_mass;
              // the residues the fragment ladder leaves for the shift, [pre_pro_index, to) in this mode's orientation: the range the DAG edge
              // tested the variable modifications on (shiftRangeEnd_ of the same two vertices and the same summed shift)
              const int len = (int)hit.getSequence().size();
              const int to = shiftRangeEnd_(t_node_spec, hi.tol_spec_map_[m], t_pro_masses, pre_node_index, pre_pro_index, node_index, pro_index, mass_shift);
              int end = m > 0 ? (to - 1) : (len - 1 - pre_pro_index);
              int start = m > 0 ? (pre_pro_index - 1) : (len - 1 - to);

              const double mod_tol = hi.tol_spec_map_[m][node_index].getIntensity();
              bool is_var = var_count != prev_var_count;
              String var_names, var_accs;
              if (is_var)
              {
                const int from = pre_pro_index;
                const int n_site = (m != 0 && pre_node_index == 0) ? pre_pro_index : -1;
                const int c_site = (m == 0 && pre_node_index == 0) ? pre_pro_index : (m == 2 && node_index == (int)t_node_spec.size() - 1) ? pro_index - 1 : -1;
                int oriented_site = -1;
                var_names = variableModificationNames_(mod_mass, mod_tol, oriented_seq, from, to, n_site, c_site, var_accs, oriented_site);
                // The edge fitted a combination on these residues; a name is missing only if the shift measured here (from the previous
                // shift as observed) is outside the tolerance that the edge's shift met (from the previous shift as searched: the nominal
                // mass of an earlier variable modification, or a mode-2 blind shift moved to the precursor mass). It is then reported as a
                // blind shift, though scored and counted as a variable one.
                if (var_names.empty())
                {
                  is_var = false;
                  if (debug) std::cout << hit.getAccession() << "\tvariable shift " << mod_mass << " without a fitting name, reported as blind" << std::endl;
                }
                else if (oriented_site >= 0)         // one possible residue: localized there
                {
                  const int pos = m > 0 ? oriented_site : len - 1 - oriented_site;
                  start = pos - 1;
                  end = pos;
                }
              }
              mod_starts.push_back(start + 1);
              mod_ends.push_back(end);
              mod_masses.push_back(mod_mass);
              mod_tols.push_back(mod_tol);
              mod_vars.push_back(is_var ? 1 : 0);
              mod_var_names.push_back(var_names);
              mod_var_accs.push_back(var_accs);
            }

            if (debug)
            {
              std::cout << hit.getAccession() << "\tmode\t" << m << "\tinput pre\t" << given_precursor_mass_ << "\tcal pre\t"
                        << std::to_string(hi.calculated_precursor_mass_) << "\tscore\t" << getScore_(*iter) << "\t ni: "
                        << node_index << "\t pi: " << pro_index
                        << "\tin\t" << t_node_spec.size() << "\t" << t_pro_masses.size() << "\tmasses\t" << t_pro_masses.back() << "\t"
                        << std::to_string(t_pro_masses[pro_index]) << "\t" << std::to_string(t_pro_masses.back() - t_pro_masses[pro_index]) << "\t"
                        << std::to_string(t_node_spec[node_index].getMZ()) << " node score " << t_node_spec[node_index].getIntensity()
                        << "\t tolspec: " << hi.tol_spec_map_[m][node_index].getIntensity() << " tol: " << t_node_spec[node_index].getMZ() / 1e5 << "\t"
                        << std::to_string(mass_shift) << "\tmod mass: " << std::to_string(mod_mass) << "\t" << mod_count << std::endl;
            }

            if (mod_count > 0 && prev_mod_count != mod_count) prev_mass_shift = mass_shift;
            prev_mod_count = mod_count;
            prev_var_count = var_count;
            pre_pro_index = pro_index;
            pre_node_index = node_index;
          }
          if (debug) { std::cout << std::endl; }
          int mode_score = getScore_(best_path[0]);
          if (m == 2 && hi.protein_end_position_ >= 0 && hi.protein_end_position_ < (int)hit.getSequence().size())
            mode_score -= truncation_penalty; // penalize C-terminal truncation in the both-termini path
          if (m == 1 && hi.protein_start_position_ >= 0 && hi.protein_end_position_ >= 0 && hi.protein_start_position_ >= hi.protein_end_position_)
          {
            // the N- and C-terminal paths overlap: keep one of them, and drop the modifications and matched positions of the other
            if (total_score > mode_score) // mode 0 wins
            {
              hi.protein_start_position_ = -1;
              erase_mods(mods_before[1], mod_starts.size());
              matched_position_map.erase(1);
              break;
            }
            else // mode 1 wins
            {
              hi.protein_end_position_ = -1;
              if (mods_before.count(0)) erase_mods(mods_before[0], mods_before[1]);
              matched_position_map.erase(0);
              total_score = mode_score;
              used_mode.pop_back();
              used_mode.push_back(m);
              break;
            }
          }
          else
          {
            total_score += mode_score;
            used_mode.push_back(m);
          }
        }

        // A mass inferred from fragments is the sum of one C- and one N-terminal fragment mass; nothing checks that both come from one
        // precursor. It is kept only with evidence in the data: the MS1 precursor mass confirms it, or the spectrum's complementary
        // fragment pairs single it out (native intact-antibody spectra: light-chain masses off by +33 to +2985 Da or by one residue).
        // Any precursor mass is also given up when the proteoform needs a blind modification beyond max_mod_mass. A rejected inferred
        // mass is replaced by the MS1 precursor mass if there is one; otherwise the terminal result is kept with an undetermined mass.
        bool reject_mass = false;
        if (hi.mode_ == 2 && total_score > 0)
        {
          // a blind modification beyond max_mod_mass is outside the search's own condition
          for (Size i = 0; i < mod_masses.size() && ! reject_mass; i++)
            reject_mass = mod_vars[i] == 0 && std::abs(mod_masses[i]) > max_mod_mass_ + std::max(0.05, i < mod_tols.size() ? mod_tols[i] : 0.0);
          const bool given_mass_agrees
            = given_precursor_mass_ > 0 && std::abs(given_precursor_mass_ - hi.calculated_precursor_mass_) <= tol_ * given_precursor_mass_;
          if (! reject_mass && precursor_by_fragment && ! given_mass_agrees) // the MS1 precursor mass confirms it
            reject_mass = ! hasComplementaryPairSupport_(sorted_masses, hi.calculated_precursor_mass_);
          // Unshifted suffix fragment masses measured from the C-terminus cannot come from a proteoform that ends before it. A precursor
          // mass of a co-isolated species makes the path end where that mass fits and absorb the rest with blind shifts: on Eclip_1480
          // scan 8324, YCII 1-76 [+28.0, -484.8] scored 37, the terminal result 1-98 with 21 y ions 228.
          const int len = (int)hit.getSequence().size();
          if (! reject_mass && terminal_end == len && terminal_c_ladder >= min_c_term_ladder && hi.protein_end_position_ < len
              && total_score < terminal_score)
            reject_mass = true;
        }
        // A mass that gives no path is given up as well when it is not the MS1 precursor mass as given (inferred from fragments, or moved by
        // an isotope offset), so that the MS1 mass is tried as it is: on Eclip_1480 scan 4032 (RL27 without the initiator Met, MS1 mass
        // 1.8 ppm from it) a fragment pair inferred the mass 1 Da lower, mode 2 found no path, and the proteoform was reported with an
        // undetermined mass.
        const bool moved = hi.mode_ == 2 && given_precursor_mass_ > 0 && ! retried_with_given_mass
                           && (precursor_by_fragment || std::abs(hi.calculated_precursor_mass_ - given_precursor_mass_) > tol_ * given_precursor_mass_);
        if (moved && total_score <= 0) reject_mass = true;
        if (reject_mass)
        {
          const bool inferred = precursor_by_fragment || moved;
          total_score = terminal_score;
          used_mode = terminal_used_mode;
          mod_starts = terminal_mod_starts;
          mod_ends = terminal_mod_ends;
          mod_masses = terminal_mod_masses;
          mod_tols = terminal_mod_tols;
          mod_vars = terminal_mod_vars;
          mod_var_names = terminal_mod_var_names;
          mod_var_accs = terminal_mod_var_accs;
          matched_position_map = terminal_matched_position_map;
          hi.protein_start_position_ = terminal_start;
          hi.protein_end_position_ = terminal_end;
          best_path_map.erase(2);
          all_path_map.erase(2);
          precursor_by_fragment = false;
          hi.calculated_precursor_mass_ = -1;
          if (inferred && ! retried_with_given_mass && given_precursor_mass_ > 0)
          {
            retried_with_given_mass = true; // mode 2 again, without inference: it takes the given precursor mass
            // as it is: the isotope offset with the most overlapping fragment masses took 1 Da off the right MS1 mass of RL27 (Eclip_1480
            // scan 4032), and the retry failed the same way
            hi.precursor_as_given_ = true;
            hi.mode_ = 1;
            continue;
          }
          break;
        }
      }

      if (terminal_kept && total_score <= 0 && terminal_score > 0)
      {
        total_score = terminal_score;
        used_mode = terminal_used_mode;
        mod_starts = terminal_mod_starts;
        mod_ends = terminal_mod_ends;
        mod_masses = terminal_mod_masses;
        mod_tols = terminal_mod_tols;
        mod_vars = terminal_mod_vars;
        mod_var_names = terminal_mod_var_names;
        mod_var_accs = terminal_mod_var_accs;
        matched_position_map = terminal_matched_position_map;
        hi.protein_start_position_ = terminal_start;
        hi.protein_end_position_ = terminal_end;
        best_path_map.erase(2);
        hi.calculated_precursor_mass_ = -1;
        precursor_by_fragment = false;
      }
      if (hi.protein_start_position_ >= 0 && hi.protein_end_position_ >= 0 && hi.protein_start_position_ >= hi.protein_end_position_) { continue; }
      if (used_mode.empty()) continue;

      const auto t_mod_masses = mod_masses, t_mod_tols = mod_tols;
      const auto t_mod_starts = mod_starts, t_mod_ends = mod_ends, t_mod_vars = mod_vars;
      const auto t_mod_var_names = mod_var_names, t_mod_var_accs = mod_var_accs;

      mod_masses.clear();
      mod_starts.clear();
      mod_ends.clear();
      mod_tols.clear();
      mod_vars.clear();
      mod_var_names.clear();
      mod_var_accs.clear();

      for (int k = 0; k < (int)t_mod_masses.size(); k++)
      {
        if (hi.protein_start_position_ >= 0 && t_mod_starts[k] < hi.protein_start_position_) continue;
        if (hi.protein_end_position_ >= 0 && t_mod_ends[k] > hi.protein_end_position_) continue;
        mod_masses.push_back(t_mod_masses[k]);
        mod_starts.push_back(t_mod_starts[k]);
        mod_ends.push_back(t_mod_ends[k]);
        mod_tols.push_back(t_mod_tols[k]);
        mod_vars.push_back(t_mod_vars[k]);
        mod_var_names.push_back(t_mod_var_names[k]);
        mod_var_accs.push_back(t_mod_var_accs[k]);
      }
      std::vector<String> mod_ids, mod_accs;

      for (int k = 0; k < (int)mod_masses.size(); k++)
      {
        if (mod_vars[k])
        {
          mod_ids.push_back(mod_var_names[k]);
          mod_accs.push_back(mod_var_accs[k]);
          continue;
        }
        auto mod_mass = mod_masses[k];
        auto iter = candidate_blind_mod_map_.lower_bound(mod_mass - mod_tols[k]);
        String mod_id = "";
        String mod_acc = "";
        std::set<int> mod_int_acc;
        while (iter != candidate_blind_mod_map_.end())
        {
          double diff = iter->first - mod_mass;
          if (diff > mod_tols[k]) break;
          if (diff > -mod_tols[k])
          {
            for (const auto& mod : iter->second)
            {
              if (mod_int_acc.find(mod.getUniModRecordId()) != mod_int_acc.end()) continue;
              mod_int_acc.insert(mod.getUniModRecordId());
              mod_acc += std::to_string(mod.getUniModRecordId()) + ",";
              mod_id += mod.getId() + ",";
            }
          }
          iter++;
        }
        mod_ids.push_back(mod_id);
        mod_accs.push_back(mod_acc);
      }
      // remove unmatched tags.
      std::set<int> to_exclude_tag_indices, matched_positions;

      for (int m = (used_mode.back() == 2 ? 2 : 0); m <= (used_mode.back() == 2 ? 2 : 1); m++)
      {
        std::vector<Size> best_path;
        const auto& t_pro_masses = hi.pro_mass_map_[m];
        if (std::find(used_mode.begin(), used_mode.end(), m) == used_mode.end() || best_path_map.empty() || best_path_map.find(m) == best_path_map.end()
            || best_path_map[m].empty())
          ;
        else
          best_path = best_path_map[m];

        if (std::find(used_mode.begin(), used_mode.end(), m) != used_mode.end())
        {
          for (int pos : matched_position_map[m])
            matched_positions.insert(pos);
        }
        for (int j = 0; j < (int)matched_tags.size(); j++) // for each tag
        {
          const auto& tag = matched_tags[j];

          if ((tag.getNtermMass() > 0 && m == 0) || (tag.getCtermMass() > 0 && m == 1)) { continue; }
          bool tag_matched = false;
          for (auto iter = best_path.rbegin(); iter != best_path.rend(); iter++) // compare against each path
          {
            auto node_index = getNodeIndex_(*iter, t_pro_masses.size());
            double node_mz = hi.node_spec_map_[m][node_index].getMZ();

            if (tag.getNtermMass() > 0)
            {
              for (const auto& shift : prefix_shifts_)
              {
                double t_mass = tag.getNtermMass() - shift;
                if (std::abs(t_mass - node_mz) > 1.1) continue;
                tag_matched = true;
                break;
              }
            }
            else
            {
              for (const auto& shift : suffix_shifts_)
              {
                double t_mass = tag.getCtermMass() - shift;
                if (m == 2 && hi.calculated_precursor_mass_ > 0) t_mass = hi.calculated_precursor_mass_ - i2f_mass - t_mass;
                if (std::abs(t_mass - node_mz) > 1.1) continue;
                tag_matched = true;
                break;
              }
            }

            if (tag_matched)
            {
              std::vector<int> positions;
              std::vector<double> masses;
              String seq = hit.getSequence();

              if (hi.protein_end_position_ >= 0) seq = seq.substr(0, hi.protein_end_position_);
              if (hi.protein_start_position_ >= 0) seq = seq.substr(hi.protein_start_position_);

              FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(positions, masses, max_mod_mass_ * max_blind_mod_cntr_ + 1, seq, tag);
              tag_matched = ! positions.empty();
              break;
            }
          }
          if (! tag_matched) to_exclude_tag_indices.insert(tag_indices[j]);
        }
      }

      std::vector<int> refined_tag_indices;
      for (auto index : tag_indices)
      {
        if (to_exclude_tag_indices.find(index) != to_exclude_tag_indices.end()) continue;
        refined_tag_indices.push_back(index);
      }

      if (refined_tag_indices.empty() || total_score <= 0) continue;

      hi.protein_start_position_ += hi.protein_start_position_ >= 0 ? 1 : 0;
      hit.setMetaValue(mi_modification_ids, mod_ids); // TODO matching masses vs. all masses?
      hit.setMetaValue(mi_modification_accs, mod_accs);
      hit.setMetaValue(mi_modifications, mod_masses);
      hit.setMetaValue(mi_modification_starts, mod_starts);
      hit.setMetaValue(mi_modification_ends, mod_ends);
      hit.setMetaValue(mi_blind_mod_count, (int)std::count(mod_vars.begin(), mod_vars.end(), 0));
      hit.setMetaValue(mi_modification_variable, mod_vars); // 1: a variable (searched) modification, 0: a blind shift
      {
        // the variable modification names of the path ("A+B," per shift, alternatives separated by ","), for FLASHTnT's penalty learning
        std::vector<String> var_entries;
        for (Size k = 0; k < mod_vars.size(); k++)
        {
          if (mod_vars[k]) var_entries.push_back(mod_var_names[k]);
        }
        hit.setMetaValue(mi_variable_modifications, var_entries);
      }
      hit.setMetaValue(mi_matched_aa, matched_positions.size());
      hit.setMetaValue(mi_tag_indices, refined_tag_indices);

      double protein_len = hit.getSequence().size();
      if (hi.protein_end_position_ > 0) { protein_len -= (protein_len - hi.protein_end_position_); }
      if (hi.protein_start_position_ > 0) { protein_len -= hi.protein_start_position_ - 1; }

      hit.setCoverage((double)matched_positions.size() / protein_len);
      hit.setScore(total_score);
      hit.setMetaValue(mi_start_position, hi.protein_start_position_);
      hit.setMetaValue(mi_end_position, hi.protein_end_position_);
      hit.setMetaValue(mi_given_mass, given_precursor_mass_);
      hit.setMetaValue(mi_mass, hi.calculated_precursor_mass_);
      hit.setMetaValue(mi_rt, dspec.getOriginalSpectrum().getRT());
      hit.setMetaValue(mi_num_mass, dspec.size());
      hit.setMetaValue(mi_precursor_score,
                       dspec.getPrecursorPeakGroup().getQscore2D()); // dspec.getPrecursorPeakGroup().getChargeSNR(dspec.getPrecursor().getCharge())
      hit.setMetaValue(mi_precursor_snr,
                       dspec.getPrecursorPeakGroup().getSNR()); // it should be charge SNR in theory. But mass SNR was written for ease of coding.
      hit.setMetaValue(mi_mass_by_fragments, precursor_by_fragment ? 1 : 0);
      // hit.setMetaValue("Proforma", string)
  #pragma omp critical
      found.emplace_back(i, hit);
    }
    endProgress();

    // Hits finish in thread order; order them by their candidate rank so that the result (and, with a single hit per spectrum, the
    // choice among equal scores) does not depend on the schedule.
    std::sort(found.begin(), found.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    for (auto& [index, hit] : found)
    {
      if (! multiple_hits_per_spec && ! proteoform_hits_.empty()) // the best score, the earliest candidate among equal ones
      {
        if (proteoform_hits_.back().getScore() >= hit.getScore()) continue;
        proteoform_hits_.pop_back();
      }
      proteoform_hits_.push_back(std::move(hit));
    }
  }

  void FLASHExtenderAlgorithm::constructDAG_(std::set<Size>& sinks,
                                             HitInformation& hi,
                                             const std::vector<std::vector<int>>& tag_edges,
                                             int max_mod_cntr_for_last_mode,
                                             bool use_tags)
  {
    Size src = getVertex_(0, 0, 0, 0, 0, hi.pro_mass_map_[hi.mode_].size());
    hi.visited_.insert(src);
    std::set<Size> visited_tag_edges;
    std::map<Size, std::set<std::pair<double, double>>> sink_map;
    std::unordered_map<Size, std::unordered_map<Size, int>> node_max_score_map; // node, cumulative mass, score

    findPathsAlongTagEndPoints(visited_tag_edges, hi, sink_map, src, 0, 0, node_max_score_map, tag_edges, max_mod_cntr_for_last_mode, use_tags);

    for (const auto& sink : sink_map)
    {
      sinks.insert(sink.first);
    }
  }

  void FLASHExtenderAlgorithm::findPathsAlongTagEndPoints(std::set<Size>& visited_tag_edges,
                                                          HitInformation& hi,
                                                          std::map<Size, std::set<std::pair<double, double>>>& sinks,
                                                          Size vertex,
                                                          double truncation_mass,
                                                          double cumulative_shift,
                                                          std::unordered_map<Size, std::unordered_map<Size, int>>& node_max_score_map,
                                                          const std::vector<std::vector<int>>& tag_edges,
                                                          int max_mod_cntr_for_last_mode,
                                                          bool use_tags)
  {
    const auto& pro_masses = *hi.cur_pro_masses_;
    int node_index = getNodeIndex_(vertex, pro_masses.size());
    int pro_index = getProIndex_(vertex, pro_masses.size());

    int tag_start_index = -1;
    int tag_end_index = -1;

    // find tag start end positions on the DAG
    const auto& tag_pro_starts = tag_edges[0];
    const auto& tag_pro_ends = tag_edges[1];
    const auto& tag_node_starts = tag_edges[2];
    const auto& tag_node_ends = tag_edges[3];

    for (int i = 0; i < (int)tag_pro_starts.size(); i++)
    {
      if (tag_start_index < 0 && tag_node_starts[i] == node_index && tag_pro_starts[i] == pro_index) { tag_start_index = i; }
      if (tag_end_index < 0 && tag_node_ends[i] == node_index && tag_pro_ends[i] == pro_index) { tag_end_index = i; }
    }

    const Size src = getVertex_(0, 0, 0, 0, 0, pro_masses.size());

    if (tag_start_index >= 0) // within tag
    {
      int node_end = -1;
      Size i = tag_start_index;
      int pro_end = -1;
      while ((i < tag_node_starts.size()) && (tag_node_starts[i] == node_index) && (tag_pro_starts[i] == pro_index))
      {
        if (pro_end < tag_pro_ends[i]) pro_end = tag_pro_ends[i];
        i++;
      }
      i = tag_start_index;
      while ((i < tag_node_starts.size()) && (tag_node_starts[i] == node_index) && (tag_pro_starts[i] == pro_index))
      {
        if (pro_end == tag_pro_ends[i]) { node_end = std::max(tag_node_ends[i], node_end); }
        i++;
      }
      std::map<Size, std::set<std::pair<double, double>>> next_vertices;

      findSubPathsBetweenTagEndPoints(next_vertices, hi, vertex, node_end, pro_end, 0, truncation_mass, cumulative_shift, node_max_score_map,
                                      max_mod_cntr_for_last_mode);

      std::map<std::set<std::pair<double, double>>, Size> mass_sink;
      for (const auto& [next_vertex, next_cumulative_shift] : next_vertices)
      {
        if (next_vertex == src) continue;
        if (mass_sink.find(next_cumulative_shift) == mass_sink.end() || getScore_(mass_sink[next_cumulative_shift]) < getScore_(next_vertex))
          mass_sink[next_cumulative_shift] = next_vertex;
      }

      for (const auto& [mass_set, next_vertex] : mass_sink)
      {
        for (const auto& masses : mass_set)
        {
          const auto& [t, c] = masses;
          findPathsAlongTagEndPoints(visited_tag_edges, hi, sinks, next_vertex, t, c, node_max_score_map, tag_edges, max_mod_cntr_for_last_mode,
                                     use_tags);
        }
      }
    }

    if (vertex == src || tag_end_index >= 0) // between tag.
    {
      std::set<Size> reachable_vertices;

      for (Size tag_index = 0; tag_index < tag_node_starts.size(); tag_index++) // for all reachable tag starting point, run extension
      {
        int node_start = tag_node_starts[tag_index];
        int pro_start = tag_pro_starts[tag_index];
        int num_blind_mod = getBlindModNumber_(vertex);
        int num_var_mod = getVarModNumber_(vertex);
        //if (hi.mode_ == 1) std::cout << pro_start << " " << max_node_start << " " << node_index   << " " << node_start << std::endl;// 0 -1
        if (pro_index > pro_start) continue; // protein end point for src no directional tag present. Maybe speed up ?
        if (node_start >= 0)
        {
          if (node_index > node_start) //
            continue;
        }
        else if (vertex != src) //
          continue;

        bool is_visited_start = false;

        for (int nbm = 0; nbm <= num_blind_mod; nbm++)
        {
          for (int nvm = 0; nvm <= num_var_mod; nvm++)
          {
            if (visited_tag_edges.find(getVertex_(node_start >= 0 ? node_start : 0, pro_start, 0, nbm, nvm, pro_masses.size()))
                != visited_tag_edges.end())
            {
              is_visited_start = true;
              break;
            }
          }
          if (is_visited_start) break;
        }

        if (is_visited_start) continue;
        visited_tag_edges.insert(getVertex_(node_start >= 0 ? node_start : 0, pro_start, 0, num_blind_mod, num_var_mod, pro_masses.size()));

        std::map<Size, std::set<std::pair<double, double>>> next_vertices;
        findSubPathsBetweenTagEndPoints(next_vertices, hi, vertex, node_start, pro_start, 0, truncation_mass, cumulative_shift, node_max_score_map,
                                        max_mod_cntr_for_last_mode);

        if (node_start < 0)
        {
          sinks = next_vertices;
          return;
        }

        std::map<std::set<std::pair<double, double>>, Size> mass_sink;

        for (const auto& [next_vertex, next_cumulative_shift] : next_vertices)
        {
          if (next_vertex == src) continue;
          if (mass_sink.find(next_cumulative_shift) == mass_sink.end() || getScore_(mass_sink[next_cumulative_shift]) < getScore_(next_vertex))
            mass_sink[next_cumulative_shift] = next_vertex;
        }

        for (const auto& [mass_set, next_vertex] : mass_sink)
        {
          for (const auto& masses : mass_set)
          {
            const auto& [t, c] = masses;
            findPathsAlongTagEndPoints(visited_tag_edges, hi, sinks, next_vertex, t, c, node_max_score_map, tag_edges, max_mod_cntr_for_last_mode,
                                       use_tags);
          }
          reachable_vertices.insert(next_vertex);
        }
      }
      if ((vertex != src || ! use_tags) && reachable_vertices.empty())
      {
        if (hi.mode_ != 2)
        {
          findSubPathsBetweenTagEndPoints(sinks, hi, vertex, -1, -1, use_tags ? 1e5 : 0, truncation_mass, cumulative_shift, node_max_score_map,
                                          max_mod_cntr_for_last_mode);
        }
        else
        {
          for (int j = 0; j < (int)pro_masses.size(); j++)
          {
            if (std::abs(hi.calculated_precursor_mass_ + truncation_mass - pro_masses[j])
                > max_mod_mass_ * (max_blind_mod_cntr_ - getBlindModNumber_(vertex)) + varModRoom_(max_var_mod_cntr_ - getVarModNumber_(vertex)))
              continue;

            findSubPathsBetweenTagEndPoints(sinks, hi, vertex, hi.cur_node_spec_->size() - 1, j, 0, truncation_mass, cumulative_shift,
                                            node_max_score_map, max_mod_cntr_for_last_mode);
          }
        }
      }
    }
  }

  void FLASHExtenderAlgorithm::findSubPathsBetweenTagEndPoints(std::map<Size, std::set<std::pair<double, double>>>& sinks,
                                                               HitInformation& hi,
                                                               Size start_vertex,
                                                               int end_node_index,
                                                               int end_pro_index,
                                                               int diagonal_counter,
                                                               double truncation_mass,
                                                               double cumulative_mod_mass,
                                                               std::unordered_map<Size, std::unordered_map<Size, int>>& node_max_score_map,
                                                               int max_blind_mod_cntr_for_last_mode)
  {
    // TODO N term mod vs. 1st amino acid mod distinction

    if (! hi.visited_.contains(start_vertex)) return;
    const auto& pro_masses = *hi.cur_pro_masses_;
    const auto pro_mass_size = pro_masses.size();
    int max_blind_mod_cntr = max_blind_mod_cntr_for_last_mode >= 0 ? max_blind_mod_cntr_for_last_mode : max_blind_mod_cntr_;
    int start_node_index = getNodeIndex_(start_vertex, pro_mass_size);
    int start_pro_index = getProIndex_(start_vertex, pro_mass_size);
    int start_score = getScore_(start_vertex);
    int start_num_blind_mod = getBlindModNumber_(start_vertex);
    int start_num_var_mod = getVarModNumber_(start_vertex);
    // no modification left: only unmodified stretches
    if (start_num_blind_mod == max_blind_mod_cntr && start_num_var_mod >= max_var_mod_cntr_) diagonal_counter = 1e5;
    const double var_mod_room = varModRoom_(max_var_mod_cntr_ - start_num_var_mod);

    const auto src = getVertex_(0, 0, 0, 0, 0, pro_mass_size);
    const auto& node_spec = *hi.cur_node_spec_;
    const auto& tol_spec = *hi.cur_tol_spec_;

    if (end_pro_index < 0) //
    {
      end_pro_index = 0;
      if (hi.protein_end_position_ >= 0)
        while (end_pro_index < pro_mass_size
               && pro_masses[end_pro_index] - pro_masses[hi.protein_end_position_ - 1] < max_mod_mass_ * (max_blind_mod_cntr - start_num_blind_mod) + var_mod_room - 1)
          end_pro_index++;
      else
        end_pro_index = ((start_num_blind_mod < max_blind_mod_cntr || start_num_var_mod < max_var_mod_cntr_) && diagonal_counter == 0)
                          ? start_pro_index + max_extension_stretch_
                          : ((int)pro_mass_size - 1); // if sink is not specified, stretch up to max_extension_stretch_ amino acids.
      end_pro_index = std::min(end_pro_index, (int)pro_mass_size - 1);
    }
    // make the range of truncation well...  make use of the positional information
    if (start_vertex == src)
    {
      double max_shift = max_mod_mass_ * (max_var_mod_cntr_ + max_blind_mod_cntr - start_num_blind_mod) + 1.1;
      double pms = hi.protein_start_position_ >= 0? pro_masses[hi.protein_start_position_] : -1;
      double delta = end_node_index >= 0? node_spec[end_node_index].getMZ() - pro_masses[end_pro_index] : 1e100;

      for (int pro_i = start_pro_index + 1; pro_i <= end_pro_index; pro_i++) // change later
      {
        double pm = pro_masses[pro_i];
        if (hi.mode_ == 0)
        {
          if(hi.c_best_shift_ - max_shift > pm) continue;
          if(hi.c_best_shift_ + max_shift < pm) break;
        }
        else if (hi.mode_ == 1)
        {
          if(hi.n_best_shift_ - max_shift > pm) continue;
          if(hi.n_best_shift_ + max_shift < pm) break;
        }

        if (pms >= 0
            && pm - pms
                 > max_shift)
          break;
        if (pms >= 0
            && pms - pm
                 > max_shift)
          continue;

        if (delta < 1e100 && std::abs(delta + pm) > max_shift)
        {
          continue;
        }

        Size vertex2 = getVertex_(0, pro_i, -truncation_penalty, 0, 0, pro_mass_size); // penalize terminal truncation
        bool connected = hi.dag_.addEdge(vertex2, start_vertex, hi.visited_);

        if (vertex2 >= hi.dag_.size() || ! connected) continue;
        findSubPathsBetweenTagEndPoints(sinks, hi, vertex2, end_node_index, end_pro_index, diagonal_counter, pro_masses[pro_i], cumulative_mod_mass,
                                        node_max_score_map, max_blind_mod_cntr_for_last_mode);
      }
    }
    // double start_delta_mass = start_node_mass - pro_masses[start_pro_index];
    if (end_node_index >= 0)
    {
      double end_node_mass = node_spec[end_node_index].getMZ();
      double end_delta_mass = end_node_mass - pro_masses[end_pro_index];
      double margin = tol_spec[end_node_index].getIntensity(); // tol_spec[end_node_index].getIntensity();
      if (std::abs(end_delta_mass - cumulative_mod_mass + truncation_mass) > max_mod_mass_ * (max_blind_mod_cntr - start_num_blind_mod) + var_mod_room + margin)
      {
        return;
      }
      if (std::abs(end_delta_mass - cumulative_mod_mass + truncation_mass) > margin)
      {
        if (diagonal_counter > 0) return; //
      }
      else
        diagonal_counter = 1e5; // if the start and end points make a diagonal line, go through the diagonal line.
    }
    bool same_score_node = false;
    Size key_to_score_map = start_pro_index * pro_mass_size + end_pro_index; // - truncation_mass; //
    for (int nbm = 0; nbm <= start_num_blind_mod; nbm++)
    {
      for (int nvm = 0; nvm <= start_num_var_mod; nvm++)
      {
        Size zero_score_vertex = getVertex_(start_node_index, start_pro_index, 0, nbm, nvm, pro_mass_size);
        auto it_vertex = node_max_score_map.find(zero_score_vertex);
        if (it_vertex != node_max_score_map.end())
        {
          auto& kts = it_vertex->second;
          auto it_score = kts.find(key_to_score_map);
          if (it_score != kts.end())
          {
            if (it_score->second > start_score) { return; }
            else if (key_to_score_map != 0 && it_score->second == start_score)
              same_score_node = true;
          }
        }
      }
    }
    if (start_node_index == end_node_index && start_pro_index == end_pro_index && start_vertex != src)
    {
      sinks[start_vertex].insert(std::pair<double, double> {truncation_mass, cumulative_mod_mass});
      return;
    }

    if (end_node_index < 0)
    {
      if (start_vertex != src) { sinks[start_vertex].insert(std::pair<double, double> {truncation_mass, cumulative_mod_mass}); }
    }
    else if (start_node_index == node_spec.size() - 1)
    {
      sinks[start_vertex].insert(std::pair<double, double> {truncation_mass, cumulative_mod_mass});
      return;
    }
    else if (start_node_index > end_node_index || start_pro_index > end_pro_index)
      return;
    if (start_num_blind_mod > 0 && same_score_node) return; // exclude truncation

    node_max_score_map[getVertex_(start_node_index, start_pro_index, 0, start_num_blind_mod, start_num_var_mod, pro_mass_size)][key_to_score_map]
      = start_score;

    for (int node_i = start_node_index + 1; node_i <= (end_node_index < 0 ? ((int)node_spec.size() - 1) : end_node_index); node_i++)
    {
      int score = start_score + (int)node_spec[node_i].getIntensity();
      double t_node_mass = node_spec[node_i].getMZ();
      double t_margin = tol_spec[node_i].getIntensity();

      const double max_mod_with_margin = std::max(max_mod_mass_, var_mod_combinations_.empty() ? 0.0 : std::max(std::abs(var_mod_combinations_.front().mass),
                                                                                                                std::abs(var_mod_combinations_.back().mass)))
                                         + t_margin;

      for (int pro_i = start_pro_index; pro_i <= end_pro_index; pro_i++)
      {
        const double pro_mass = pro_masses[pro_i]; // 캐싱
        const double t_delta_mass = t_node_mass - pro_mass;
        const double delta_delta = t_delta_mass - cumulative_mod_mass + truncation_mass;

        if (delta_delta > max_mod_with_margin) continue;
        if (delta_delta < -max_mod_with_margin) break;

        int num_var_mod = start_num_var_mod;
        int num_blind_mod = start_num_blind_mod;
        int next_score = score;
        const int series = (*hi.cur_node_series_)[node_i];
        if (series >= 0) next_score += hi.cur_bond_weight_[series][pro_i]; // cleavage propensity of the bond this mass falls on
        double next_cumulative_mod_mass = cumulative_mod_mass;

        const double abs_delta_delta = std::abs(delta_delta); //

        if (abs_delta_delta > t_margin)
        {
          if (pro_i == start_pro_index) continue;
          if (abs_delta_delta < 0.036386 - t_margin) continue;
          if (abs_delta_delta > 0.036386 + t_margin && abs_delta_delta < 0.947630 - t_margin) continue;

          // variable modifications: a known shift on the residues of this stretch (start_pro_index, pro_i] that can carry it, counted apart
          // from blind ones. Those residues end before a mass in between that already carries the whole shift (shiftRangeEnd_); the report
          // of the path names the modification on the same residues. A shift that no combination fits there is a blind one, with the blind
          // penalty and count: tested on the whole stretch, such shifts were scored and counted as variable and then reported as blind, up
          // to 3 blind shifts with max_mod_count 2 (Eclip_1480 scans 8498 RL11, 6613 CSPA).
          bool var_search_success = false;
          if (! var_mod_combinations_.empty() && num_var_mod < max_var_mod_cntr_ && diagonal_counter == 0)
          {
            // terminal residues this stretch reaches: the path start (N-terminal for prefix paths, C-terminal for suffix paths) and, in
            // mode 2, the precursor mass node at the C-terminus
            const int n_site = (hi.mode_ != 0 && start_node_index == 0) ? start_pro_index : -1;
            const int c_site = (hi.mode_ == 0 && start_node_index == 0) ? start_pro_index
                               : (hi.mode_ == 2 && node_i == (int)node_spec.size() - 1) ? pro_i - 1 : -1;
            int range_end = -1; // computed once, for the first combination that fits the whole stretch
            auto iter = std::lower_bound(var_mod_combinations_.begin(), var_mod_combinations_.end(), delta_delta - t_margin,
                                         [](const VariableModificationCombination& c, double m) { return c.mass < m; });
            for (; iter != var_mod_combinations_.end() && iter->mass <= delta_delta + t_margin; ++iter)
            {
              if (num_var_mod + (int)iter->mods.size() > max_var_mod_cntr_) continue;
              // the whole stretch first (cheap; a combination that does not fit it fits no part of it), then the residues that can carry it
              if (! combinationFits_(*iter, hi.cur_oriented_seq_, start_pro_index, pro_i, n_site, c_site)) continue;
              if (range_end < 0) range_end = shiftRangeEnd_(node_spec, tol_spec, pro_masses, start_node_index, start_pro_index, node_i, pro_i, t_delta_mass);
              if (range_end < pro_i && ! combinationFits_(*iter, hi.cur_oriented_seq_, start_pro_index, range_end, n_site, c_site)) continue;
              var_search_success = true;
              num_var_mod += (int)iter->mods.size();
              next_cumulative_mod_mass = cumulative_mod_mass + iter->mass;
              for (int mod : iter->mods) next_score -= var_mod_penalties_[mod];
              break;
            }
          }

          if (! var_search_success)
          {
            num_blind_mod++;
            if (num_blind_mod > max_blind_mod_cntr) continue;
            if (abs_delta_delta > max_mod_mass_ + t_margin) continue; // variable combinations may reach further; blind shifts may not
            if (diagonal_counter > 0) continue;

            next_cumulative_mod_mass = t_delta_mass + truncation_mass;

            if (hi.mode_ == 2)
            {
              int start_pi = std::max(pro_i, hi.protein_end_position_ - 5);
              int end_pi = std::min(hi.protein_end_position_ + 5, (int)pro_mass_size);
              for (int pi = start_pi; pi < end_pi; pi++)
              {
                const double corrected_mod_mass = hi.calculated_precursor_mass_ - i2f_mass - pro_masses[pi] + truncation_mass;
                const double delta = corrected_mod_mass - next_cumulative_mod_mass;
                if (delta > 1.1) continue;
                if (delta < -1.1) break;
                next_cumulative_mod_mass = corrected_mod_mass;
                break;
              }
              if (std::abs(t_delta_mass - next_cumulative_mod_mass + truncation_mass) > t_margin) continue;
            }

            next_score -= 1 + 2 * (multi_ion_score + FLASHTaggerAlgorithm::max_node_score);
            auto iter = candidate_blind_mod_map_.lower_bound(delta_delta - t_margin);
            if (iter == candidate_blind_mod_map_.end() || std::abs(delta_delta - iter->first) > t_margin)
              next_score -= multi_ion_score + FLASHTaggerAlgorithm::max_node_score;
          }
        }

        next_score = std::min(next_score, max_path_score_);
        if (next_score < min_path_score_) continue;

        const Size next_vertex = getVertex_(node_i, pro_i, next_score, num_blind_mod, num_var_mod, pro_mass_size);
        if (next_vertex >= hi.dag_.size()) continue;
        if (next_score >= max_path_score_ && !sinks.empty() && hi.visited_.contains(next_vertex)) continue;
        if (!hi.dag_.addEdge(next_vertex, start_vertex, hi.visited_)) continue;

        int next_diagonal_counter = diagonal_counter;
        if (diagonal_counter > 0) next_diagonal_counter--;
        else if (num_blind_mod != start_num_blind_mod || num_var_mod != start_num_var_mod)
          next_diagonal_counter = 1;

        findSubPathsBetweenTagEndPoints(
          sinks, hi, next_vertex, end_node_index, end_pro_index, next_diagonal_counter,
          truncation_mass, next_cumulative_mod_mass, node_max_score_map, max_blind_mod_cntr_for_last_mode
        );
      }
    }
  }
} // namespace OpenMS
