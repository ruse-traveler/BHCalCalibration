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
#include <TCut.h>
#include <TROOT.h>
#include <TSystem.h>
// calibration class
#include "src/BHCalCalibration.hh"

// load libraries
R__LOAD_LIBRARY(./src/BHCalCalibration_cc.so)

using namespace std;



void DoBHCalCalibration() {

  // io parameters
  const pair<string, string> sOutput = {
    "test.root",
    "ntCalibratedEnergies"
  };
  const pair<string, string> sInput  = {
    "../performance/eicrecon_output/single_particles/merged/forPerformanceStudy.withIndividualECalLayers_includedEPar7.e110th45n20Kneu.d20m7y2023.plugin.root",
    "JCalibrateHCalWithImaging/ntForCalibration"
  };

  // how to use each leaf
  enum Usage { Tuple, Train, Target, Watcher };

  // tuple parameters
  vector<pair<int, string>> vecTrainAndTupleVars = {
    {Usage::Target, "ePar"},
    {Usage::Tuple,  "fracParVsLeadBHCal"},
    {Usage::Tuple,  "fracParVsLeadBEMC"},
    {Usage::Tuple,  "fracParVsSumBHCal"},
    {Usage::Tuple,  "fracParVsSumBEMC"},
    {Usage::Tuple,  "fracLeadBHCalVsBEMC"},
    {Usage::Tuple,  "fracSumBHCalVsBEMC"},
    {Usage::Train,  "eLeadBHCal"},
    {Usage::Train,  "eLeadBEMC"},
    {Usage::Tuple,  "eSumBHCal"},
    {Usage::Tuple,  "eSumBEMC"},
    {Usage::Tuple,  "diffLeadBHCal"},
    {Usage::Tuple,  "diffLeadBEMC"},
    {Usage::Tuple,  "diffSumBHCal"},
    {Usage::Tuple,  "diffSumBEMC"},
    {Usage::Train,  "nHitsLeadBHCal"},
    {Usage::Train,  "nHitsLeadBEMC"},
    {Usage::Tuple,  "nClustBHCal"},
    {Usage::Tuple,  "nClustBEMC"},
    {Usage::Tuple,  "hLeadBHCal"},
    {Usage::Tuple,  "hLeadBEMC"},
    {Usage::Tuple,  "fLeadBHCal"},
    {Usage::Tuple,  "fLeadBEMC"},
    {Usage::Tuple,  "eLeadImage"},
    {Usage::Train,  "eSumImage"},
    {Usage::Tuple,  "eLeadSciFi"},
    {Usage::Train,  "eSumSciFi"},
    {Usage::Tuple,  "nClustImage"},
    {Usage::Tuple,  "nClustSciFi"},
    {Usage::Tuple,  "hLeadImage"},
    {Usage::Tuple,  "hLeadSciFi"},
    {Usage::Tuple,  "fLeadImage"},
    {Usage::Tuple,  "fLeadSciFi"},
    {Usage::Train,  "eSumSciFiLayer1"},
    {Usage::Train,  "eSumSciFiLayer2"},
    {Usage::Train,  "eSumSciFiLayer3"},
    {Usage::Train,  "eSumSciFiLayer4"},
    {Usage::Train,  "eSumSciFiLayer5"},
    {Usage::Train,  "eSumSciFiLayer6"},
    {Usage::Train,  "eSumSciFiLayer7"},
    {Usage::Train,  "eSumSciFiLayer8"},
    {Usage::Train,  "eSumSciFiLayer9"},
    {Usage::Train,  "eSumSciFiLayer10"},
    {Usage::Train,  "eSumSciFiLayer11"},
    {Usage::Train,  "eSumSciFiLayer12"},
    {Usage::Train,  "eSumImageLayer1"},
    {Usage::Train,  "eSumImageLayer2"},
    {Usage::Train,  "eSumImageLayer3"},
    {Usage::Train,  "eSumImageLayer4"},
    {Usage::Train,  "eSumImageLayer5"},
    {Usage::Train,  "eSumImageLayer6"}
  };

  // tmva parameters
  const string sLoader  = "TMVADir";
  const string sFactory = "TMVARegression";

  // general tmva options
  const string sFactOpt  = "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression";
  const string sTrainOpt = "nTrain_Regression=1000:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V";
  const string sReadOpt  = "!Color:!Silent";
  const TCut   cSelect   = "eSumBHCal>0";
  const bool   addSpecs  = false;

  // tmva methods and method-specific options
  const vector<pair<string, string>> vecMethodsAndOpts = {
    {"LD",   ""},
    {"MLP",  ""},
    {"BDTG", ""}
  };

  // sort tuple leaves
  vector<string> vecTupleLeaves;
  vector<string> vecTmvaTrainers;
  vector<string> vecTmvaTargets;
  vector<string> vecTmvaSpectators;
  for (const auto trainAndTupleVars : vecTrainAndTupleVars) {
    switch (trainAndTupleVars.first) {
      case Usage::Train:
        vecTmvaTrainers.push_back(trainAndTupleVars.second);
        break;
      case Usage::Target:
        vecTmvaTargets.push_back(trainAndTupleVars.second);
        break;
      case Usage::Watcher:
        vecTmvaSpectators.push_back(trainAndTupleVars.second);
        break;
    }
    vecTupleLeaves.push_back(trainAndTupleVars.second);
  }

  // run calibration workflow
  BHCalCalibration* calibrator = new BHCalCalibration(sFactory, sLoader, sOutput.first);
  calibrator -> SetInput(sInput.first, sInput.second);
  calibrator -> SetTupleLeaves(vecTupleLeaves);
  calibrator -> SetTmvaOpts(sFactOpt, sTrainOpt, sReadOpt, addSpecs);
  calibrator -> SetTmvaArgs(vecTmvaTrainers, vecTmvaTargets, vecTmvaSpectators, cSelect);
  calibrator -> SetTmvaMethods(vecMethodsAndOpts);
  calibrator -> Init();
  calibrator -> Train();
  calibrator -> Apply();
  calibrator -> End();
  return;

}

// end ------------------------------------------------------------------------
