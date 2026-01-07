from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from libcpp cimport bool

from String cimport *
from FASTAFile cimport *
from Param cimport *
from FLASHHelperClasses cimport *

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>" namespace "OpenMS":

    cdef cppclass FLASHIda:
        # wrap-doc:
        #   FLASHIda class for real-time deconvolution and proteoform identification.
        #
        #   This class contains functions to perform deconvolution for spectra received
        #   from mass spectrometers, as well as proteoform identification using the
        #   FLASHTnT workflow.

        # Constructor taking a string argument
        FLASHIda(char* arg) except + nogil

        # Copy constructor
        FLASHIda(FLASHIda &) except + nogil

        # Python-friendly method for getting top fragment matches
        int getTopFragmentMatchesPy(String & protein_sequence,
                                    int n,
                                    libcpp_vector[double] & masses,
                                    libcpp_vector[double] & qscores,
                                    libcpp_vector[int] & charges,
                                    libcpp_vector[double] & window_starts,
                                    libcpp_vector[double] & window_ends,
                                    libcpp_vector[int] & is_b_ions,
                                    libcpp_vector[int] & fragment_indices) except + nogil
            # wrap-doc:
            #   Get top fragment ions matching the protein sequence.
            #
            #   Returns top n fragment ions sorted by qscore.
            #   Requires deconvolveMS2() to be called first.
            #
            #   :param protein_sequence: the protein sequence to match against
            #   :param n: maximum number of ions to return
            #   :param masses: output monoisotopic masses
            #   :param qscores: output quality scores
            #   :param charges: output representative charges
            #   :param window_starts: output isolation window start m/z
            #   :param window_ends: output isolation window end m/z
            #   :param is_b_ions: output 1 if b-ion, 0 if y-ion
            #   :param fragment_indices: output 1-based fragment index (b3=3, y5=5)
            #   :returns: number of matches found

        # Python-friendly method for getting PTM ambiguity enclosing ions
        int getAmbiguityEnclosingIonsPy(String & protein_sequence,
                                        int n,
                                        libcpp_vector[double] & masses,
                                        libcpp_vector[double] & qscores,
                                        libcpp_vector[int] & charges,
                                        libcpp_vector[double] & window_starts,
                                        libcpp_vector[double] & window_ends,
                                        libcpp_vector[int] & is_b_ions,
                                        libcpp_vector[int] & fragment_indices) except + nogil
            # wrap-doc:
            #   Get unique fragment ions that enclose PTM ambiguity regions.
            #
            #   Identifies PTM sites and returns the best fragment ions that bracket
            #   ambiguity regions. Returns deduplicated ions sorted by qscore.
            #   Useful for MS3 targeting to resolve PTM localization.
            #   Requires deconvolveMS2() to be called first.
            #
            #   :param protein_sequence: the protein sequence to analyze
            #   :param n: maximum number of ions to return
            #   :param masses: output monoisotopic masses
            #   :param qscores: output quality scores
            #   :param charges: output representative charges
            #   :param window_starts: output isolation window start m/z
            #   :param window_ends: output isolation window end m/z
            #   :param is_b_ions: output 1 if b-ion, 0 if y-ion
            #   :param fragment_indices: output 1-based fragment index (b3=3, y5=5)
            #   :returns: number of enclosing ions found

        # Python-friendly method for getting terminal fragment ions
        int getTerminalFragmentIonsPy(String & protein_sequence,
                                      int n,
                                      libcpp_vector[double] & masses,
                                      libcpp_vector[double] & qscores,
                                      libcpp_vector[int] & charges,
                                      libcpp_vector[double] & window_starts,
                                      libcpp_vector[double] & window_ends,
                                      libcpp_vector[int] & is_b_ions,
                                      libcpp_vector[int] & fragment_indices) except + nogil
            # wrap-doc:
            #   Get terminal (innermost) fragment ions.
            #
            #   Returns fragment ions that extend furthest toward the opposite terminus.
            #   Output is interleaved: [top_b, top_y, 2nd_b, 2nd_y, ...]
            #   Requires deconvolveMS2() to be called first.
            #
            #   :param protein_sequence: the protein sequence to match against
            #   :param n: maximum number of ions to return
            #   :param masses: output monoisotopic masses
            #   :param qscores: output quality scores
            #   :param charges: output representative charges
            #   :param window_starts: output isolation window start m/z
            #   :param window_ends: output isolation window end m/z
            #   :param is_b_ions: output 1 if b-ion, 0 if y-ion
            #   :param fragment_indices: output 1-based fragment index (b3=3, y5=5)
            #   :returns: number of ions found

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>" namespace "OpenMS::FLASHIda":

    cdef cppclass FLASHIdaTagMatch "OpenMS::FLASHIda::TagMatch":
        # wrap-doc:
        #   Structure representing a tag match to a protein database entry.

        FLASHIdaTagMatch() except + nogil
        FLASHIdaTagMatch(FLASHIdaTagMatch &) except + nogil

        String tag_sequence
        double n_term_mass
        double c_term_mass
        double tag_score
        int protein_index
        String protein_accession
        int match_position
        double flanking_mass_diff
