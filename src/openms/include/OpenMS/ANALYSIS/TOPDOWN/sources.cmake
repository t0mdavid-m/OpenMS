### the directory name
set(directory include/OpenMS/ANALYSIS/TOPDOWN)

### list all header files of the directory here
set(sources_list_h
        DeconvolvedSpectrum.h
        SpectralDeconvolution.h
        FLASHDeconvAlgorithm.h
        FLASHHelperClasses.h
        FLASHExtenderAlgorithm.h
        FLASHIda.h
        FLASHIda/Config.h
        FLASHIda/Deconvolution.h
        FLASHIda/Exploration.h
        FLASHIda/FAIMS.h
        FLASHIda/FragmentAnalysis.h
        FLASHIda/Logger.h
        FLASHIda/MS3FragmentMatcher.h
        FLASHIda/PrecursorSelection.h
        FLASHIda/Quantification.h
        FLASHIda/ScanCommand.h
        FLASHIda/ScanCommandQueue.h
        FLASHIdaBridgeFunctions.h
        FLASHGappedTaggerAlgorithm.h
        FLASHTaggerAlgorithm.h
        FLASHTnTAlgorithm.h
        MassFeatureTrace.h
        PeakGroup.h
        PeakGroupScoring.h
        Qvalue.h
        TopDownIsobaricQuantification.h
)

### add path to the filenames
set(sources_h)
foreach (i ${sources_list_h})
    list(APPEND sources_h ${directory}/${i})
endforeach (i)

### source group definition
source_group("Header Files\\OpenMS\\ANALYSIS\\TOPDOWN" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

