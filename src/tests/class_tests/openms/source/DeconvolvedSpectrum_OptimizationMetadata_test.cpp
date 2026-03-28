// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/OptimizationMetadata.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

using namespace OpenMS;

START_TEST(DeconvolvedSpectrum_OptimizationMetadata, "$Id$")

/////////////////////////////////////////////////////////////

// P2-U01: Default-constructed DeconvolvedSpectrum has no metadata
START_SECTION(hasOptimizationMetadata_default_false)
{
  DeconvolvedSpectrum ds(1); // scan_number = 1
  TEST_EQUAL(ds.hasOptimizationMetadata(), false)
}
END_SECTION

// P2-U02: getOrCreateOptimizationMetadata creates metadata
START_SECTION(getOrCreateOptimizationMetadata_creates_and_returns_true)
{
  DeconvolvedSpectrum ds(1);
  TEST_EQUAL(ds.hasOptimizationMetadata(), false)
  OptimizationMetadata& meta = ds.getOrCreateOptimizationMetadata();
  TEST_EQUAL(ds.hasOptimizationMetadata(), true)
  TEST_NOT_EQUAL(ds.getOptimizationMetadata(), (const OptimizationMetadata*)nullptr)
  // Suppress unused variable warning
  (void)meta;
}
END_SECTION

// P2-U03: All metadata fields carry correct defaults
START_SECTION(metadata_field_defaults)
{
  DeconvolvedSpectrum ds(1);
  OptimizationMetadata& meta = ds.getOrCreateOptimizationMetadata();
  TEST_EQUAL(meta.group_id, 0)
  TEST_EQUAL(meta.variant_index, -1)
  TEST_EQUAL(meta.total_variants, 0)
  TEST_EQUAL(meta.is_best_variant, false)
  TEST_EQUAL(meta.rank, 0)
  TEST_EQUAL(meta.msn_level_optimized, 0)
  TEST_EQUAL(meta.parent_tracking_id, 0)
  TEST_REAL_SIMILAR(meta.collision_energy, 0.0)
  TEST_REAL_SIMILAR(meta.isolation_width, 0.0)
  TEST_STRING_EQUAL(meta.activation_type, "")
  TEST_REAL_SIMILAR(meta.precursor_mass, 0.0)
  TEST_EQUAL(meta.precursor_charge, 0)
  TEST_REAL_SIMILAR(meta.fragmentation_quality_score, -1.0)
  TEST_EQUAL(meta.fragment_count, 0)
  TEST_EQUAL(meta.exploration_scans, 0)
}
END_SECTION

// P2-U04: toSpectrum serializes metadata via setMetaValue when present
START_SECTION(toSpectrum_serializes_metadata_via_setMetaValue)
{
  DeconvolvedSpectrum ds(2); // scan_number = 2

  // toSpectrum() accesses peak_groups_[0], so we need at least one PeakGroup
  PeakGroup pg;
  ds.push_back(pg);

  OptimizationMetadata& meta = ds.getOrCreateOptimizationMetadata();
  meta.group_id = 42;
  meta.collision_energy = 25.0;
  meta.is_best_variant = true;
  meta.fragmentation_quality_score = 0.87;
  meta.precursor_mass = 15432.5;

  MSSpectrum out_spec = ds.toSpectrum(1);

  TEST_EQUAL((int)out_spec.getMetaValue("optimization_group_id"), 42)
  TEST_REAL_SIMILAR((double)out_spec.getMetaValue("optimization_collision_energy"), 25.0)
  TEST_STRING_EQUAL((std::string)out_spec.getMetaValue("optimization_is_best_variant"), "true")
  TEST_REAL_SIMILAR((double)out_spec.getMetaValue("optimization_quality_score"), 0.87)
  TEST_REAL_SIMILAR((double)out_spec.getMetaValue("optimization_precursor_mass"), 15432.5)
}
END_SECTION

// P2-U05: toSpectrum without metadata sets no optimization metavalues
START_SECTION(toSpectrum_without_metadata_sets_no_optimization_metavalues)
{
  DeconvolvedSpectrum ds(1); // scan_number = 1
  TEST_EQUAL(ds.hasOptimizationMetadata(), false)

  // toSpectrum() accesses peak_groups_[0], so we need at least one PeakGroup
  PeakGroup pg;
  ds.push_back(pg);

  MSSpectrum out_spec = ds.toSpectrum(1);

  TEST_EQUAL(out_spec.metaValueExists("optimization_group_id"), false)
  TEST_EQUAL(out_spec.metaValueExists("optimization_collision_energy"), false)
  TEST_EQUAL(out_spec.metaValueExists("optimization_is_best_variant"), false)
  TEST_EQUAL(out_spec.metaValueExists("optimization_quality_score"), false)
  TEST_EQUAL(out_spec.metaValueExists("optimization_precursor_mass"), false)
}
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
