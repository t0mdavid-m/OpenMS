// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

// pyOpenMS nanobind bindings
// Domain: flashida
//
// Exposes FLASHTnTAlgorithm and PeakGroup for the FLASHIda RL training pipeline.

#include "all_casters.h"
#include "binding_utils.h"

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ProteinHit.h>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_pyopenms_flashida, m) {
    m.doc() = "pyOpenMS FLASHIda bindings (FLASHTnTAlgorithm, PeakGroup)";

    // -----------------------------------------------------------------------
    // PeakGroup::TargetDecoyType enum
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::PeakGroup::TargetDecoyType>(m, "PeakGroupTargetDecoyType",
        "Specifies if a PeakGroup is a target (0), noise decoy (1), or signal decoy (2)")
        .value("target", OpenMS::PeakGroup::TargetDecoyType::target)
        .value("noise_decoy", OpenMS::PeakGroup::TargetDecoyType::noise_decoy)
        .value("signal_decoy", OpenMS::PeakGroup::TargetDecoyType::signal_decoy)
        .export_values();

    // -----------------------------------------------------------------------
    // PeakGroup
    // Accessors needed by the RL training pipeline.
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakGroup>(m, "PeakGroup",
        R"doc(
Deconvolved peak group (mass + charge envelope + isotope pattern).

A mass contains multiple LogMz peaks of different charges and isotope indices.
PeakGroup is the set of such peaks representing a single monoisotopic mass and
contains the quality features used for PeakGroupScoring.
DeconvolvedSpectrum consists of PeakGroups.
)doc")
        .def(nb::init<>())
        .def(nb::init<int, int, bool>(), "min_abs_charge"_a, "max_abs_charge"_a, "is_positive"_a,
             "Constructor specifying charge range and ionisation mode")
        .def(nb::init<const OpenMS::PeakGroup&>())
        .def("__copy__", [](const OpenMS::PeakGroup& self) { return OpenMS::PeakGroup(self); })
        .def("__deepcopy__", [](const OpenMS::PeakGroup& self, nb::dict) {
            return OpenMS::PeakGroup(self);
        }, "memo"_a)

        // --- basic mass / intensity ---
        .def("getMonoMass", &OpenMS::PeakGroup::getMonoMass,
             "Returns the monoisotopic mass")
        .def("getIntensity", &OpenMS::PeakGroup::getIntensity,
             "Returns the summed intensity over all signal peaks")

        // --- per-charge queries ---
        .def("getChargeIntensity", &OpenMS::PeakGroup::getChargeIntensity, "abs_charge"_a,
             "Returns the summed intensity for the given absolute charge state")
        .def("getChargeSNR", &OpenMS::PeakGroup::getChargeSNR, "abs_charge"_a,
             "Returns the SNR for the given absolute charge state")
        .def("getChargeIsotopeCosine", &OpenMS::PeakGroup::getChargeIsotopeCosine, "abs_charge"_a,
             "Returns the isotope cosine score for the given absolute charge state")

        // --- global quality scores ---
        .def("getIsotopeCosine", &OpenMS::PeakGroup::getIsotopeCosine,
             "Returns the cosine score between averagine and the observed isotope pattern")
        .def("getSNR", &OpenMS::PeakGroup::getSNR,
             "Returns the total SNR")
        .def("getQscore", &OpenMS::PeakGroup::getQscore,
             "Returns the quality score")
        .def("getChargeScore", &OpenMS::PeakGroup::getChargeScore,
             "Returns the charge fit score")
        .def("getAvgPPMError", &OpenMS::PeakGroup::getAvgPPMError,
             "Returns the average mass ppm error")

        // --- charge range / representative charge ---
        .def("getRepAbsCharge", &OpenMS::PeakGroup::getRepAbsCharge,
             "Returns the representative charge (the one maximising SNR)")
        .def("getAbsChargeRange", [](const OpenMS::PeakGroup& self) {
            auto [mn, mx] = self.getAbsChargeRange();
            return nb::make_tuple(mn, mx);
        }, "Returns the absolute charge range as a (min_charge, max_charge) tuple")

        // --- targeting ---
        .def("isTargeted", &OpenMS::PeakGroup::isTargeted,
             "Returns True if this peak group has been targeted")
        .def("getTargetDecoyType", &OpenMS::PeakGroup::getTargetDecoyType,
             "Returns the TargetDecoyType (target / noise_decoy / signal_decoy)")

        // --- scan ---
        .def("getScanNumber", &OpenMS::PeakGroup::getScanNumber,
             "Returns the scan number this peak group was found in")
        ;

    // -----------------------------------------------------------------------
    // FLASHHelperClasses::Tag
    // Used as output of FLASHTnTAlgorithm::getTags() and input to
    // getProteoformHitsMatchedBy().
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FLASHHelperClasses::Tag>(m, "Tag",
        R"doc(
Sequence tag used by FLASHTaggerAlgorithm / FLASHTnTAlgorithm.
No mass gap is allowed within the sequence; mass-gap-containing tags must be
enumerated into multiple Tag instances externally.
)doc")
        .def(nb::init<const OpenMS::FLASHHelperClasses::Tag&>())
        .def("__copy__", [](const OpenMS::FLASHHelperClasses::Tag& self) {
            return OpenMS::FLASHHelperClasses::Tag(self);
        })
        .def("__deepcopy__", [](const OpenMS::FLASHHelperClasses::Tag& self, nb::dict) {
            return OpenMS::FLASHHelperClasses::Tag(self);
        }, "memo"_a)
        .def("getSequence", &OpenMS::FLASHHelperClasses::Tag::getSequence,
             "Returns the tag sequence (with modification annotations)")
        .def("getUppercaseSequence", &OpenMS::FLASHHelperClasses::Tag::getUppercaseSequence,
             "Returns the uppercase sequence (no modification annotations)")
        .def("getNtermMass", &OpenMS::FLASHHelperClasses::Tag::getNtermMass,
             "Returns the N-terminal flanking mass")
        .def("getCtermMass", &OpenMS::FLASHHelperClasses::Tag::getCtermMass,
             "Returns the C-terminal flanking mass")
        .def("getLength", &OpenMS::FLASHHelperClasses::Tag::getLength,
             "Returns the tag length (number of residues)")
        .def("getScore",
             nb::overload_cast<>(&OpenMS::FLASHHelperClasses::Tag::getScore, nb::const_),
             "Returns the overall tag score")
        .def("getScore",
             nb::overload_cast<int>(&OpenMS::FLASHHelperClasses::Tag::getScore, nb::const_),
             "pos"_a, "Returns the score at position pos")
        .def("getScan", &OpenMS::FLASHHelperClasses::Tag::getScan,
             "Returns the scan number from which the tag was generated")
        .def("getIndex", &OpenMS::FLASHHelperClasses::Tag::getIndex,
             "Returns the tag index")
        .def("setIndex", &OpenMS::FLASHHelperClasses::Tag::setIndex, "i"_a,
             "Sets the tag index")
        .def("getRetentionTime", &OpenMS::FLASHHelperClasses::Tag::getRetentionTime,
             "Returns the retention time of the source spectrum")
        .def("setRetentionTime", &OpenMS::FLASHHelperClasses::Tag::setRetentionTime, "rt"_a,
             "Sets the retention time")
        .def("getMzs", &OpenMS::FLASHHelperClasses::Tag::getMzs,
             "Returns the m/z values of the peaks used to generate the tag")
        .def("toString", &OpenMS::FLASHHelperClasses::Tag::toString,
             "Returns a human-readable string representation")
        .def("__repr__", [](const OpenMS::FLASHHelperClasses::Tag& self) {
            return "Tag(seq='" + std::string(self.getSequence()) + "', scan=" + std::to_string(self.getScan()) + ")";
        })
        ;

    // -----------------------------------------------------------------------
    // FLASHTnTAlgorithm
    // -----------------------------------------------------------------------
    auto flashtnt_class =
        nb::class_<OpenMS::FLASHTnTAlgorithm, OpenMS::DefaultParamHandler>(m, "FLASHTnTAlgorithm",
        R"doc(
FLASHTnT top-down identification algorithm.

FLASHTnT (FLASH Tag-and-Test) identifies proteoforms in deconvolved
MSExperiments produced by FLASHDeconv. It builds sequence tags from MS2
spectra, searches them against a FASTA proteome, and returns ProteinHit
objects representing matched proteoforms.

DefaultParamHandler
ProgressLogger
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FLASHTnTAlgorithm&>())
        .def("__copy__", [](const OpenMS::FLASHTnTAlgorithm& self) {
            return OpenMS::FLASHTnTAlgorithm(self);
        })
        .def("__deepcopy__", [](const OpenMS::FLASHTnTAlgorithm& self, nb::dict) {
            return OpenMS::FLASHTnTAlgorithm(self);
        }, "memo"_a)

        .def("run",
            [](OpenMS::FLASHTnTAlgorithm& self,
               const OpenMS::MSExperiment& exp,
               const std::vector<OpenMS::FASTAFile::FASTAEntry>& fasta) {
                nb::gil_scoped_release release;
                self.run(exp, fasta);
            },
            "exp"_a, "fasta"_a,
            R"doc(
Run identification on a deconvolved MSExperiment against a FASTA proteome.

Decoy or MS level 1 spectra are filtered out internally. Overlapping
PeakGroups in merged spectra are also removed.

:param exp:   MSExperiment produced by FLASHDeconvAlgorithm.run()
:param fasta: List of FASTAEntry objects to search against
)doc")

        .def("getProteoforms",
            [](const OpenMS::FLASHTnTAlgorithm& self) {
                std::vector<OpenMS::ProteinHit> hits;
                self.getProteoforms(hits);
                return hits;
            },
            "Returns the list of identified proteoform ProteinHit objects")

        .def("getTags",
            [](const OpenMS::FLASHTnTAlgorithm& self) {
                std::vector<OpenMS::FLASHHelperClasses::Tag> tags;
                self.getTags(tags);
                return tags;
            },
            "Returns all sequence tags generated during the run")

        .def("getProteoformHitsMatchedBy",
            [](const OpenMS::FLASHTnTAlgorithm& self,
               const OpenMS::FLASHHelperClasses::Tag& tag) {
                std::vector<OpenMS::ProteinHit> hits;
                self.getProteoformHitsMatchedBy(tag, hits);
                return hits;
            },
            "tag"_a,
            "Returns the proteoform ProteinHits that are matched/supported by the given Tag")

        .def("getMaxTotalModificationMass",
            &OpenMS::FLASHTnTAlgorithm::getMaxTotalModificationMass,
            "Returns the maximum total modification mass allowed during search")
        ;

    // FLASHTnTAlgorithm inherits both DefaultParamHandler and ProgressLogger.
    // DefaultParamHandler is declared as nanobind base above; ProgressLogger
    // methods are injected via the template helper to avoid multiple-inheritance issues.
    def_ProgressLogger<OpenMS::FLASHTnTAlgorithm>(flashtnt_class);

} // NB_MODULE
