// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jaekwan Kim$
// $Authors: Jaekwan Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/FORMAT/MzMLFile.h>
//#include <OpenMS/FORMAT/MzMLFile.h>
///////////////////////////


using namespace OpenMS;
using namespace std;

START_TEST(FLASHTaggerAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FLASHTaggerAlgorithm* ptr = 0;
FLASHTaggerAlgorithm* null_ptr = 0;
START_SECTION(FLASHTaggerAlgorithm())
{
  ptr = new FLASHTaggerAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FLASHTaggerAlgorithm())
{
  delete ptr;
}
END_SECTION

PeakMap map;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SIGMAmAb_ECDHCD_nativeTD_1_deconv.mzML"), map);

std::vector<DeconvolvedSpectrum> deconvolved_spectra;
double tol(10);
// run FLASHDeconvAlgorithm here!


// collect statistics for information
for (int index = 0; index < map.size(); index++)
{

  auto spec = map[index];
  if (spec.getMSLevel() != 2) continue;
  int scan = FLASHDeconvAlgorithm::getScanNumber(map, index);
  DeconvolvedSpectrum dspec(scan);
  dspec.setOriginalSpectrum(spec);
  String deconv_meta_str = spec.getMetaValue("DeconvMassInfo").toString();
  int tol_loc_s = deconv_meta_str.find("tol=") + 4;
  int tol_loc_e = deconv_meta_str.find(";", tol_loc_s);
  tol = stod(deconv_meta_str.substr(tol_loc_s, tol_loc_e - tol_loc_s));

  int q_loc_s = deconv_meta_str.find("qscore=") + 7;
  int q_loc_e = deconv_meta_str.find(";", q_loc_s);
  auto q_str = deconv_meta_str.substr(q_loc_s, q_loc_e - q_loc_s);
  Size pos = 0;
  std::vector<double> qscores;
  while (true)
  {
    Size pos_t = q_str.find(",", pos);
    if (pos_t == String::npos) break;
    auto token = q_str.substr(pos, pos_t - pos);
    qscores.push_back(stod(token));
    pos = pos_t + 1;
  }

  int s_loc_s = deconv_meta_str.find("snr=") + 4;
  int s_loc_e = deconv_meta_str.find(";", s_loc_s);
  auto s_str = deconv_meta_str.substr(s_loc_s, s_loc_e - s_loc_s);
  pos = 0;
  std::vector<float> snrs;
  while (true)
  {
    Size pos_t = s_str.find(",", pos);
    if (pos_t == String::npos) break;
    auto token = s_str.substr(pos, pos_t - pos);
    snrs.push_back(stof(token));
    pos = pos_t + 1;
  }

  for (int i = 0; i < spec.size(); i++)
  {
    PeakGroup peak;
    peak.setQscore(qscores[i]);
    peak.setSNR(snrs[i]);
    peak.setMonoisotopicMass(spec[i].getMZ());
    peak.setScanNumber(scan);
    dspec.push_back(peak);
  }
  dspec.sort();
  deconvolved_spectra.push_back(dspec);
}
// Run tagger
FLASHTaggerAlgorithm tagger;




START_SECTION(run())
{
  tagger.run(deconvolved_spectra, tol);
  TEST_EQUAL(tagger.getTags().size(),264)
}
END_SECTION


START_SECTION(runMatching())
{
  tagger.runMatching(OPENMS_GET_TEST_DATA_PATH("uniprot-mg1655-filtered-reviewed_yes.fasta"));
  TEST_EQUAL(tagger.getProteinHits().size(), 2)
}
END_SECTION


START_SECTION(getProteinHits())
{
  std::string sequence1 = "EVQLV";
  std::string sequence2 = "GPSLS";
  TEST_EQUAL(tagger.getProteinHits()[0].getSequence().substr(0, 5), sequence1)
  TEST_EQUAL(tagger.getProteinHits()[1].getSequence().substr(0, 5), sequence2)
}
END_SECTION


START_SECTION(getProteinHits(const FLASHDeconvHelperStructs::Tag& tag)) 
{
  TEST_EQUAL(tagger.getProteinHits(tagger.getTags()[0]).size(), 1)
  TEST_EQUAL(tagger.getProteinHits(tagger.getTags()[1]).size(), 0)
}
END_SECTION

START_SECTION(getProteinHits(const FLASHDeconvHelperStructs::Tag& tag)) 
{
  TEST_EQUAL(tagger.getProteinHits(tagger.getTags()[0]).size(), 1)
  TEST_EQUAL(tagger.getProteinHits(tagger.getTags()[1]).size(), 0)
}
END_SECTION


START_SECTION(getTags(const ProteinHit& hit))
{
  std::string sequence3 = "GQGTMVTVSS";
  std::string sequence4 = "SVTVMTGQGW";
  TEST_EQUAL(tagger.getTags(tagger.getProteinHits()[0])[0].getSequence(), sequence3)
  TEST_EQUAL(tagger.getTags(tagger.getProteinHits()[1])[1].getSequence(), sequence4)
}
END_SECTION

START_SECTION(getTags())
{
  TEST_EQUAL(tagger.getTags()[0].getScore(), 18) 
  TEST_EQUAL(tagger.getTags()[1].getScore(), 18)
}
END_SECTION

START_SECTION(getProteinIndex())
 {
 
   TEST_EQUAL(tagger.getProteinIndex(tagger.getProteinHits()[0]), 0) 
   TEST_EQUAL(tagger.getProteinIndex(tagger.getProteinHits()[1]), 1)
 }
END_SECTION

START_SECTION(getTagIndex())
 {
 
   TEST_EQUAL(tagger.getTagIndex(tagger.getTags()[0]), 0) 
   TEST_EQUAL(tagger.getTagIndex(tagger.getTags()[1]), 1)
   TEST_EQUAL(tagger.getTagIndex(tagger.getTags()[2]), 2) 
   TEST_EQUAL(tagger.getTagIndex(tagger.getTags()[3]), 3)
   TEST_EQUAL(tagger.getTagIndex(tagger.getTags()[4]), 4) 


 }
END_SECTION

START_SECTION(getMatchedPositions()) 
{

    TEST_EQUAL(tagger.getMatchedPositions(tagger.getProteinHits()[0], tagger.getTags()[0])[0], 111) 

  } END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST