from libcpp.vector cimport vector as libcpp_vector
from MSExperiment cimport *
from FASTAFile cimport *
from ProteinHit cimport ProteinHit
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from FLASHHelperClasses cimport Tag


cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>" namespace "OpenMS":

    cdef cppclass FLASHTnTAlgorithm(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #   DefaultParamHandler
        #   ProgressLogger

        # default constructor
        FLASHTnTAlgorithm() except + nogil
        # copy constructor
        FLASHTnTAlgorithm(FLASHTnTAlgorithm &) except + nogil

        void run(MSExperiment & map, libcpp_vector[FASTAEntry] & fasta_entry) except + nogil

        void getProteoformHitsMatchedBy(Tag & tag, libcpp_vector[ProteinHit] & hits) except + nogil

        void getTags(libcpp_vector[Tag] & tags) except + nogil

        void getProteoforms(libcpp_vector[ProteinHit] & hits) except + nogil

        double getMaxTotalModificationMass() except + nogil
