# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Kyowon Jeong $
# $Authors: Kyowon Jeong $
# --------------------------------------------------------------------------

# CRoaring 5.2.2 (https://github.com/RoaringBitmap/CRoaring), the amalgamated release files roaring.c, roaring.h and roaring.hh:
# compressed bitmaps of 32-bit integers, used for the protein window index of FLASHGappedTaggerAlgorithm.
add_library(CRoaring STATIC roaring.c)
target_include_directories(CRoaring SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(CRoaring PROPERTIES POSITION_INDEPENDENT_CODE ON C_STANDARD 11 C_STANDARD_REQUIRED ON)

get_filename_component(CROARING_FOLDER "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(APPEND ${TPLICENSEFILE} "\nCRoaring:\n")
cat(${CROARING_FOLDER}/LICENSE ${TPLICENSEFILE})
