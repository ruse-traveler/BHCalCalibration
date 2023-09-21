// ----------------------------------------------------------------------------
// 'DoBHCalCalibration.cc'
// Derek Anderson
// 09.14.2023
//
// Driver macro to run ePIC BHCal calibration workflow.
// ----------------------------------------------------------------------------

// c utilities
#include <string>
#include <vector>
#include <utility>
// root utilites
#include <TROOT.h>
#include <TSystem.h>
// calibration class
#include "src/BHCalCalibration.hh"

// load libraries
R__LOAD_LIBRARY(./src/BHCalCalibration_cc.so)

using namespace std;



void DoBHCalCalibration() {

  // load library
  //gROOT -> ProcessLine(".L ./src/BHCalCalibration.cc++");
  //gSystem -> Load("./src/BHCalCalibration_cc.so");

  // io parameters
  const string           sOutput = "test.root";
  const array<string, 2> sInput  = {
    "../performance/eicrecon_output/single_particles/merged/forPerformanceStudy.withIndividualECalLayers_includedEPar7.e110th45n20Kneu.d20m7y2023.plugin.root",
    "JCalibrateWithImaging/ntToCalibrate"
  };

  // tuple parameters
  vector<string> vecTupleVars = {
    "ePar",
    "fracParVsLeadBHCal",
    "fracParVsLeadBEMC",
    "fracParVsSumBHCal",
    "fracParVsSumBEMC",
    "fracLeadBHCalVsBEMC",
    "fracSumBHCalVsBEMC",
    "eLeadBHCal",
    "eLeadBEMC",
    "eSumBHCal",
    "eSumBEMC",
    "diffLeadBHCal",
    "diffLeadBEMC",
    "diffSumBHCal",
    "diffSumBEMC",
    "nHitsLeadBHCal",
    "nHitsLeadBEMC",
    "nClustBHCal",
    "nClustBEMC",
    "hLeadBHCal",
    "hLeadBEMC",
    "fLeadBHCal",
    "fLeadBEMC",
    "eLeadImage",
    "eSumImage",
    "eLeadSciFi",
    "eSumSciFi",
    "nClustImage",
    "nClustSciFi",
    "hLeadImage",
    "hLeadSciFi",
    "fLeadImage",
    "fLeadSciFi",
    "eSumSciFiLayer1",
    "eSumSciFiLayer2",
    "eSumSciFiLayer3",
    "eSumSciFiLayer4",
    "eSumSciFiLayer5",
    "eSumSciFiLayer6",
    "eSumSciFiLayer7",
    "eSumSciFiLayer8",
    "eSumSciFiLayer9",
    "eSumSciFiLayer10",
    "eSumSciFiLayer11",
    "eSumSciFiLayer12",
    "eSumImageLayer1",
    "eSumImageLayer2",
    "eSumImageLayer3",
    "eSumImageLayer4",
    "eSumImageLayer5",
    "eSumImageLayer6"
  };

  // tmva parameters
  const string sLoader("TMVADir");
  const string sFactory("TMVARegression");

  // run calibration workflow
  BHCalCalibration* calibrator = new BHCalCalibration(sFactory, sLoader, sOutput);
  calibrator -> SetInput(sInput[0], sInput[1]);
  calibrator -> SetTupleArgs(vecTupleVars);
  calibrator -> Init();
  calibrator -> Train();
  calibrator -> Apply();
  calibrator -> End();
  return;

}

// end ------------------------------------------------------------------------
