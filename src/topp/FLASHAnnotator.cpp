// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
@page TOPP_FLASHAnnotator FLASHAnnotator

  @brief FLASHAnnotator annotates the PrSMs reported by FLASHTnT on the raw MS2 spectra.

  For each PrSM of a FLASHTnT PrSM file (-out_prsm of FLASHTnT), the theoretical isotope envelope of every terminal
  fragment of the reported proteoform is placed at every plausible charge on the centroided raw MS2 spectrum (profile
  spectra are centroided with PeakPickerHiRes), the way ExD Viewer annotates, and the envelopes that match are kept
  (see FLASHAnnotatorAlgorithm). FLASHAnnotator reports the bond coverage and a PrSM score that is calibrated against
  the chance match rate of mass-shifted decoy fragments in the same spectrum, so that it does not grow with the
  proteoform length.

  The PrSM rows are written unchanged, followed by the columns RawScore (-log10 of the binomial tail probability of the
  matched fragment hypotheses), RawLLR (log-likelihood ratio), RawMatchedFragments, RawFragmentHypotheses, RawChanceRate,
  RawCoveredBonds, RawCoverage(%) and RawCalibrationPPM. -1 marks values that could not be computed (scan not found,
  unknown residue; for RawScore and RawChanceRate also no decoy shifts or no fragment hypothesis); RawLLR and
  RawCalibrationPPM, which can be negative, are NaN then. RawFragmentHypotheses = -1 marks a PrSM that was not annotated.
  -out_annotation lists the accepted envelopes.

  RawScore and RawLLR are evidence scores of a PrSM in its own spectrum, not identification p-values: FLASHTnT selects
  its candidates by tags and blind modifications fitted to the spectrum, so even a wrong candidate matches more fragments
  than a random proteoform. For identification, use target-decoy competition (e.g. the scores of the PrSMs of decoy
  proteins, written by FLASHTnT with -keep_decoy).

  Fragment charges are capped by the precursor charge: the charge of the raw precursor if the precursor mass of the PrSM
  (PrecursorMass) puts an isotope peak in the isolation window at that charge; otherwise the charge implied by the mass
  and the isolation target m/z if that charge does, else the larger of the two. Without a raw charge, the implied charge
  is used; without a precursor mass either, the fragment charges are capped only by algorithm:min_mass_per_charge. The
  proteoform mass is not used: it is not the precursor mass when the precursor is a complex.

  The input mzML must be the raw file FLASHDeconv was run on (FLASHTnT scan numbers refer to it). Only the MS2 spectra of
  the PrSM scans are kept in memory.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FLASHAnnotator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FLASHAnnotator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

/// keeps, of the spectra read by MzMLFile::transform, the first MS2 (or higher level) spectrum of every needed scan number
class NeededSpectraConsumer : public Interfaces::IMSDataConsumer
{
public:
  /// with @p scan_from_index, scan numbers are spectrum index + 1 unless the native ID carries one (all spectra must be read)
  NeededSpectraConsumer(const std::set<int>& needed_scans, bool scan_from_index): needed_scans_(needed_scans), scan_from_index_(scan_from_index)
  {
  }

  void setExperimentalSettings(const ExperimentalSettings& settings) override
  {
    // as FLASHDeconvAlgorithm::getScanNumber
    type_accession_ = "MS:1000768";
    if (! settings.getSourceFiles().empty() && ! settings.getSourceFiles()[0].getNativeIDTypeAccession().empty())
    {
      type_accession_ = settings.getSourceFiles()[0].getNativeIDTypeAccession();
    }
  }

  void setExpectedSize(Size, Size) override
  {
  }

  void consumeSpectrum(SpectrumType& spec) override
  {
    const int index = read_count_++;
    std::vector<String> native_ids;
    spec.getNativeID().split(",", native_ids);
    int scan = native_ids.empty() ? -1 : SpectrumLookup::extractScanNumber(native_ids.back(), type_accession_);
    if (scan < 0)
    {
      without_scan_number_ = true;
      if (! scan_from_index_) return;
      scan = index + 1;
    }
    if (spec.getMSLevel() < 2 || needed_scans_.find(scan) == needed_scans_.end() || scan_to_slot.find(scan) != scan_to_slot.end()) return;
    scan_to_slot[scan] = (int)spectra.size();
    spectra.push_back(std::move(spec));
  }

  void consumeChromatogram(ChromatogramType&) override
  {
  }

  /// the kept spectra and the scan number -> index into spectra map
  std::vector<MSSpectrum> spectra;
  std::map<int, int> scan_to_slot;

  int readCount() const
  {
    return read_count_;
  }

  /// true if a spectrum had no scan number in its native ID
  bool withoutScanNumber() const
  {
    return without_scan_number_;
  }

private:
  const std::set<int>& needed_scans_;
  bool scan_from_index_;
  String type_accession_ = "MS:1000768";
  int read_count_ = 0;
  bool without_scan_number_ = false;
};

class TOPPFLASHAnnotator : public TOPPBase
{
public:
  TOPPFLASHAnnotator():
      TOPPBase("FLASHAnnotator",
               "Annotates FLASHTnT PrSMs on the raw MS2 spectra with fragment isotope envelopes and scores them against the chance "
               "match rate.",
               false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input raw mzML file (profile or centroid) the PrSMs were identified from.", true);
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("in_prsm", "<file>", "", "Input PrSM tsv file (-out_prsm of FLASHTnT).", true);
    setValidFormats_("in_prsm", ListUtils::create<String>("tsv"));

    registerOutputFile_("out", "<file>", "", "Output PrSM tsv file: the input rows followed by the raw annotation columns.", true);
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_annotation", "<file>", "", "Output tsv file with one row per accepted fragment isotope envelope.", false);
    setValidFormats_("out_annotation", ListUtils::create<String>("tsv"));

    registerStringList_("ion_type", "<types>", ListUtils::toStringList(std::vector<std::string> {"b", "y"}),
                        "Ion types to annotate (as for FLASHTnT).", false);
    setValidStrings_("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

    registerIntOption_("precursor_charge", "<charge>", 0,
                       "Precursor charge used for all PrSMs (caps the fragment charges). 0: the charge of the raw precursor if "
                       "the PrecursorMass of the PrSM has an isotope peak in the isolation window at it; else "
                       "round(PrecursorMass / (precursor m/z - proton)) if that charge does, or if the "
                       "raw charge is not annotated; else the larger of the two; without a mass or precursor m/z, no cap other "
                       "than algorithm:min_mass_per_charge.",
                       false);
    setMinInt_("precursor_charge", 0);

    registerSubsection_("algorithm", "Annotation algorithm parameters");
  }

  Param getSubsectionDefaults_(const String& prefix) const override
  {
    if (prefix == "algorithm")
    {
      auto param = FLASHAnnotatorAlgorithm().getDefaults();
      param.remove("ion_type");
      return param;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", prefix);
  }

  static double toDouble_(String s, double default_value)
  {
    try
    {
      s.trim();
      return s.empty() ? default_value : s.toDouble();
    }
    catch (Exception::BaseException&)
    {
      return default_value;
    }
  }

  static int toInt_(String s, int default_value)
  {
    try
    {
      s.trim();
      return s.empty() ? default_value : (int)std::lround(s.toDouble());
    }
    catch (Exception::BaseException&)
    {
      return default_value;
    }
  }

  /// ';'-separated values, empty ones skipped
  static std::vector<String> splitList_(const String& s)
  {
    std::vector<String> tokens, values;
    s.split(";", tokens);
    for (auto& t : tokens)
    {
      t.trim();
      if (! t.empty()) values.push_back(t);
    }
    return values;
  }

  static double seconds_(const std::chrono::steady_clock::time_point& since)
  {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - since).count();
  }

  /// @p value as an ostream writes it, or "NaN"
  static std::string format_(double value)
  {
    if (std::isnan(value)) return "NaN";
    std::ostringstream ss;
    ss << value;
    return ss.str();
  }

  /// the appended PrSM columns and the envelope header
  static constexpr const char* prsm_columns_ = "\tRawScore\tRawLLR\tRawMatchedFragments\tRawFragmentHypotheses\tRawChanceRate\tRawCoveredBonds"
                                                "\tRawCoverage(%)\tRawCalibrationPPM";
  static constexpr const char* envelope_header_ = "PrSMIndex\tScan\tIonType\tLength\tBond\tCharge\tMonoMass\tMostAbundantMz\tMatchedPeaks"
                                                   "\tPredictedPeaks\tCosine\tPearson\tIntensity\tPPMError";

  ExitCodes main_(int, const char**) override
  {
    const auto start_time = std::chrono::steady_clock::now();
    String in_file = getStringOption_("in");
    String in_prsm_file = getStringOption_("in_prsm");
    String out_file = getStringOption_("out");
    String out_annotation_file = getStringOption_("out_annotation");
    const int fixed_precursor_charge = getIntOption_("precursor_charge");

    Param algorithm_param = getParam_().copy("algorithm:", true);
    algorithm_param.setValue("ion_type", getParam_().getValue("ion_type"));
    FLASHAnnotatorAlgorithm annotator;
    annotator.setParameters(algorithm_param);

    //-------------------------------------------------------------
    // PrSMs: columns are found by their header names
    //-------------------------------------------------------------
    auto time = std::chrono::steady_clock::now();
    std::ifstream prsm_stream(in_prsm_file);
    if (! prsm_stream) { throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, in_prsm_file); }
    String header;
    std::vector<String> lines;
    std::vector<std::vector<String>> rows;
    {
      std::string line;
      bool first = true;
      while (std::getline(prsm_stream, line))
      {
        if (! line.empty() && line.back() == '\r') line.pop_back();
        if (first)
        {
          header = line;
          first = false;
          continue;
        }
        if (String(line).trim().empty()) continue;
        lines.emplace_back(line);
        std::vector<String> fields;
        lines.back().split('\t', fields);
        rows.push_back(fields);
      }
    }
    std::vector<String> columns;
    header.split('\t', columns);
    auto column = [&columns](const String& name) -> int {
      for (Size i = 0; i < columns.size(); i++)
      {
        if (String(columns[i]).trim() == name) return (int)i;
      }
      return -1;
    };
    int col_index = column("PrSMIndex");
    if (col_index < 0) col_index = 0;
    const int col_scan = column("Scan");
    const int col_sequence = column("ProteinSequence");
    const int col_start = column("StartPosition");
    const int col_end = column("EndPosition");
    const int col_precursor_mass = column("PrecursorMass");
    const int col_mod_mass = column("ModMass");
    const int col_mod_start = column("ModStart");
    const int col_mod_end = column("ModEnd");
    if (col_scan < 0 || col_sequence < 0)
    {
      OPENMS_LOG_ERROR << "The PrSM file has no Scan or ProteinSequence column: " << in_prsm_file << endl;
      return INCOMPATIBLE_INPUT_DATA;
    }
    auto field = [](const std::vector<String>& row, int col) -> String { return (col >= 0 && col < (int)row.size()) ? row[col] : String(); };

    std::set<int> needed_scans;
    for (const auto& row : rows)
    {
      needed_scans.insert(toInt_(field(row, col_scan), -1));
    }
    const double read_time = seconds_(time);

    if (rows.empty())
    {
      // nothing to annotate: the raw file is not read
      std::ofstream out(out_file);
      if (! out) { throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_file); }
      out << header << prsm_columns_ << "\n";
      if (! out_annotation_file.empty())
      {
        std::ofstream out_annotation(out_annotation_file);
        if (! out_annotation) { throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_annotation_file); }
        out_annotation << envelope_header_ << "\n";
      }
      OPENMS_LOG_INFO << "PrSMs: 0" << endl;
      return EXECUTION_OK;
    }

    //-------------------------------------------------------------
    // raw spectra: only the MS2 spectra of the PrSM scans are kept
    //-------------------------------------------------------------
    time = std::chrono::steady_clock::now();
    MzMLFile mzml;
    mzml.setLogType(log_type_);
    PeakFileOptions options = mzml.getOptions();
    options.setMSLevels({2, 3, 4, 5});
    // chromatograms are not skipped: with setSkipChromatograms(true), MzMLHandler ends a progress it did not start
    mzml.setOptions(options);
    auto consumer = std::make_unique<NeededSpectraConsumer>(needed_scans, false);
    mzml.transform(in_file, consumer.get(), true);
    if (consumer->withoutScanNumber())
    {
      // scan numbers then come from the spectrum indices, which require all spectra
      OPENMS_LOG_INFO << "Native IDs without scan numbers: reading all spectra" << endl;
      options.setMSLevels({});
      mzml.setOptions(options);
      consumer = std::make_unique<NeededSpectraConsumer>(needed_scans, true);
      mzml.transform(in_file, consumer.get(), true);
    }
    std::vector<MSSpectrum>& raw_spectra = consumer->spectra;
    const std::map<int, int>& scan_to_slot = consumer->scan_to_slot;
    const double load_time = seconds_(time);

    //-------------------------------------------------------------
    // the spectra of the PrSMs: centroided (if profile) once per spectrum, in parallel; the raw peaks are then released
    //-------------------------------------------------------------
    time = std::chrono::steady_clock::now();
    std::vector<FLASHAnnotatorAlgorithm::PeakArrays> spectra(raw_spectra.size());
    PeakPickerHiRes picker;
    int profile_count = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : profile_count)
    for (int s = 0; s < (int)raw_spectra.size(); s++)
    {
      MSSpectrum& spec = raw_spectra[s];
      if (spec.getType(true) == SpectrumSettings::PROFILE)
      {
        MSSpectrum picked;
        picker.pick(spec, picked);
        FLASHAnnotatorAlgorithm::toPeakArrays(picked, spectra[s]);
        profile_count++;
      }
      else { FLASHAnnotatorAlgorithm::toPeakArrays(spec, spectra[s]); }
      spec.clear(false);
      spec.shrink_to_fit();
    }

    // PrSM inputs
    std::vector<FLASHAnnotatorAlgorithm::PrSM> prsms(rows.size());
    std::vector<double> prsm_masses(rows.size(), -1.0);
    double max_prsm_mass = 0;
    int missing_scan = 0;
    for (Size r = 0; r < rows.size(); r++)
    {
      const auto& row = rows[r];
      auto& prsm = prsms[r];
      prsm.sequence = field(row, col_sequence).trim();
      const int start = toInt_(field(row, col_start), 1);
      const int end = toInt_(field(row, col_end), 1);
      // FLASHTnT writes -1 for a start or end it could not determine: those fragments have no defined mass
      prsm.known_n_term = start > 0;
      prsm.known_c_term = end > 0;
      // modification positions are 1-based protein positions
      const int start_in_seq = start < 0 ? 0 : start - 1;
      const auto mod_masses = splitList_(field(row, col_mod_mass));
      const auto mod_starts = splitList_(field(row, col_mod_start));
      const auto mod_ends = splitList_(field(row, col_mod_end));
      for (Size m = 0; m < mod_masses.size() && m < mod_starts.size() && m < mod_ends.size(); m++)
      {
        FLASHAnnotatorAlgorithm::Modification mod;
        mod.mass = toDouble_(mod_masses[m], 0);
        mod.start = toInt_(mod_starts[m], 1) - 1 - start_in_seq;
        mod.end = toInt_(mod_ends[m], 1) - 1 - start_in_seq;
        prsm.modifications.push_back(mod);
      }

      // the precursor mass only: the proteoform mass is not the precursor mass when the precursor is a complex or a
      // co-isolated species, and a charge derived from it would cap the fragment charges too low
      const double mass = toDouble_(field(row, col_precursor_mass), -1);
      prsm_masses[r] = mass;
      max_prsm_mass = std::max({max_prsm_mass, mass, FLASHAnnotatorAlgorithm::getProteoformMass(prsm.sequence, prsm.modifications)});

      const auto it = scan_to_slot.find(toInt_(field(row, col_scan), -1));
      if (it == scan_to_slot.end())
      {
        missing_scan++;
        continue;
      }
      prsm.spectrum_index = it->second;
    }

    // precursor charges. The isotope table is computed once, beyond the heaviest PrSM (the annotation then reuses it).
    int charge_from_raw = 0, charge_from_prsm = 0, charge_larger = 0, charge_from_mass = 0, charge_uncapped = 0;
    std::vector<double> iso_mzs, iso_intensities;
    if (fixed_precursor_charge == 0) annotator.getIsotopeEnvelope(max_prsm_mass + 200.0, 1, iso_mzs, iso_intensities);
    // true if an isotope peak of mass at charge z lies in [low, high]
    auto in_window = [&](double mass, int z, double low, double high) -> bool {
      annotator.getIsotopeEnvelope(mass, z, iso_mzs, iso_intensities);
      return ! iso_mzs.empty() && iso_mzs.front() <= high && iso_mzs.back() >= low;
    };
    for (Size r = 0; r < rows.size(); r++)
    {
      auto& prsm = prsms[r];
      if (prsm.spectrum_index < 0) continue;
      if (fixed_precursor_charge > 0)
      {
        prsm.precursor_charge = fixed_precursor_charge;
        continue;
      }
      const auto& spec = raw_spectra[prsm.spectrum_index];
      const double mass = prsm_masses[r];
      int raw_charge = 0;
      double precursor_mz = 0, low = 0, high = 0;
      if (! spec.getPrecursors().empty())
      {
        const auto& precursor = spec.getPrecursors()[0];
        raw_charge = precursor.getCharge();
        precursor_mz = precursor.getMZ();
        low = precursor_mz - precursor.getIsolationWindowLowerOffset();
        high = precursor_mz + precursor.getIsolationWindowUpperOffset();
      }
      // the charge at which the precursor mass has the isolation target m/z
      const int implied_charge = (mass > 0 && precursor_mz > Constants::PROTON_MASS_U) ?
                                   std::max(1, (int)std::round(mass / (precursor_mz - Constants::PROTON_MASS_U))) : 0;
      if (raw_charge > 0)
      {
        // a harmonic or co-isolated raw charge does not fit the proteoform: fragments above it would never be tried
        if (implied_charge <= 0 || implied_charge == raw_charge || in_window(mass, raw_charge, low, high))
        {
          prsm.precursor_charge = raw_charge;
          charge_from_raw++;
        }
        else if (in_window(mass, implied_charge, low, high))
        {
          prsm.precursor_charge = implied_charge;
          charge_from_prsm++;
        }
        else
        {
          prsm.precursor_charge = std::max(raw_charge, implied_charge);
          charge_larger++;
        }
      }
      else if (implied_charge > 0)
      {
        prsm.precursor_charge = implied_charge;
        charge_from_mass++;
      }
      else
      {
        // no precursor charge or mass: the fragment charges are capped by their own masses only
        prsm.precursor_charge = 0;
        charge_uncapped++;
      }
    }
    const double centroid_time = seconds_(time);

    //-------------------------------------------------------------
    // annotation
    //-------------------------------------------------------------
    time = std::chrono::steady_clock::now();
    std::vector<FLASHAnnotatorAlgorithm::Result> results;
    annotator.annotate(spectra, prsms, results);
    const double annotation_time = seconds_(time);

    //-------------------------------------------------------------
    // output, in the input order
    //-------------------------------------------------------------
    time = std::chrono::steady_clock::now();
    {
      std::ofstream out(out_file);
      if (! out) { throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_file); }
      out << header << prsm_columns_ << "\n";
      for (Size r = 0; r < rows.size(); r++)
      {
        const auto& res = results[r];
        out << lines[r] << "\t" << res.score << "\t" << format_(res.llr) << "\t" << res.matched_fragments << "\t" << res.fragment_hypotheses
            << "\t" << res.chance_rate << "\t" << res.covered_bonds << "\t" << res.coverage << "\t" << format_(res.calibration_ppm) << "\n";
      }
    }
    if (! out_annotation_file.empty())
    {
      std::ofstream out(out_annotation_file);
      if (! out) { throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_annotation_file); }
      out << envelope_header_ << "\n";
      const auto& ion_types = annotator.getIonTypes();
      for (Size r = 0; r < rows.size(); r++)
      {
        const String prsm_index = field(rows[r], col_index);
        const String scan = field(rows[r], col_scan);
        std::ostringstream ss;
        for (const auto& e : results[r].envelopes)
        {
          ss << prsm_index << "\t" << scan << "\t" << ion_types[e.ion_type] << "\t" << e.length << "\t" << e.bond << "\t" << e.charge << "\t"
             << std::fixed << std::setprecision(5) << e.mono_mass << "\t" << e.most_abundant_mz << "\t" << e.matched_peaks << "\t"
             << e.predicted_peaks << "\t" << std::setprecision(4) << e.cosine << "\t" << e.pearson << "\t" << std::setprecision(1) << e.intensity
             << "\t" << std::setprecision(2) << e.ppm_error << "\n";
        }
        out << ss.str();
      }
    }
    const double write_time = seconds_(time);

    Size annotated = 0, envelopes = 0;
    for (const auto& res : results)
    {
      annotated += res.annotated ? 1 : 0;
      envelopes += res.envelopes.size();
    }
    OPENMS_LOG_INFO << "PrSMs: " << rows.size() << " (" << annotated << " annotated, " << missing_scan << " without their scan in the mzML, "
                    << (rows.size() - annotated - missing_scan) << " with an unknown residue); accepted envelopes: " << envelopes << endl;
    OPENMS_LOG_INFO << "Spectra: " << consumer->readCount() << " read, " << spectra.size() << " kept and annotated (" << profile_count
                    << " centroided)" << endl;
    if (fixed_precursor_charge == 0)
    {
      OPENMS_LOG_INFO << "Precursor charge: " << charge_from_raw << " from the raw precursor, " << charge_from_prsm
                      << " from the precursor mass and the isolation window (raw charge does not fit), " << charge_larger
                      << " the larger of both (neither fits), " << charge_from_mass << " from the precursor mass (no raw charge), "
                      << charge_uncapped << " without a precursor cap" << endl;
    }
    OPENMS_LOG_INFO << std::fixed << std::setprecision(2) << "Timing (s): PrSM file " << read_time << ", mzML load " << load_time
                    << ", centroiding " << centroid_time << ", annotation " << annotation_time << ", writing " << write_time
                    << ", total " << seconds_(start_time) << endl;
    return EXECUTION_OK;
  }
};

/// @endcond

// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  TOPPFLASHAnnotator tool;
  return tool.main(argc, argv);
}
