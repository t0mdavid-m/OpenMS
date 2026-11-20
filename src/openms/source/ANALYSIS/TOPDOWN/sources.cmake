### the directory name
set(directory source/ANALYSIS/TOPDOWN)

### list all filenames of the directory here
set(sources_list
        DeconvolvedSpectrum.cpp
        FLASHAnnotatorAlgorithm.cpp
        SpectralDeconvolution.cpp
        FLASHDeconvAlgorithm.cpp
        FLASHHelperClasses.cpp
        FLASHExtenderAlgorithm.cpp
        FLASHIda.cpp
        FLASHIda/Config.cpp
        FLASHIda/ConsoleView.cpp
        FLASHIda/Deconvolution.cpp
        FLASHIda/Exploration.cpp
        FLASHIda/EngineClock.cpp
        FLASHIda/EventLog.cpp
        FLASHIda/FAIMS.cpp
        FLASHIda/AnnotatorMetrics.cpp
        FLASHIda/BuildInfo.cpp
        FLASHIda/ChanceMatch.cpp
        FLASHIda/ChargeCheck.cpp
        FLASHIda/KineticModel.cpp
        FLASHIda/Refine.cpp
        FLASHIda/FragmentAnalysis.cpp
        FLASHIda/IdaLogger.cpp
        FLASHIda/LogFormat.cpp
        FLASHIda/MassOffset.cpp
        FLASHIda/MS3FragmentMatcher.cpp
        FLASHIda/Ms1Buffer.cpp
        FLASHIda/Ms3WitnessGate.cpp
        FLASHIda/Ms3WitnessModel.cpp
        FLASHIda/Ms3YieldModel.cpp
        FLASHIda/PrecursorBoard.cpp
        FLASHIda/PrecursorSelection.cpp
        FLASHIda/ProteoformTracker.cpp
        FLASHIda/ProteomeSearch.cpp
        FLASHIda/Quantification.cpp
        FLASHIda/ScanCommandQueue.cpp
        FLASHIda/SearchPool.cpp
        FLASHIdaBridgeFunctions.cpp
        FLASHGappedTaggerAlgorithm.cpp
        FLASHTaggerAlgorithm.cpp
        FLASHTnTAlgorithm.cpp
        FLASHTnTCandidateRanker.cpp
        FLASHQuantAlgorithm.cpp
        FLASHQuantHelper.cpp
        MassFeatureTrace.cpp
        PeakGroup.cpp
        PeakGroupScoring.cpp
        Qvalue.cpp
        TopDownIsobaricQuantification.cpp
        MassLadderAlgorithm.cpp
        ProteoformDatabase.cpp
        FixedModifications.cpp
)

### add path to the filenames
set(sources)
foreach (i ${sources_list})
    list(APPEND sources ${directory}/${i})
endforeach (i)

### Engine build 1: the CI run that built the engine (events.log [ENGINE]), on BuildInfo.cpp alone so a new run
### recompiles one small file. A digits-only id; outside CI the file reports "local".
if(DEFINED ENV{GITHUB_RUN_ID} AND "$ENV{GITHUB_RUN_ID}" MATCHES "^[0-9]+$")
  set_source_files_properties(${directory}/FLASHIda/BuildInfo.cpp PROPERTIES COMPILE_DEFINITIONS "FLASHIDA_CI_RUN=$ENV{GITHUB_RUN_ID}")
endif()

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\TOPDOWN" FILES ${sources})

