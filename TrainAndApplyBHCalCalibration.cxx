// ----------------------------------------------------------------------------
// 'TrainAndApplyBHCalCalibration.C'
// Derek Anderson
// 09.13.2023
//
// Macro to test ePIC BHCal calibration workflow.  Ingests
// TNtuple summarizes info from BHCal and BECal and trains/
// applies TMVA model based on specified parameters.
// ----------------------------------------------------------------------------

// c utilities
#include <map>
#include <array>
#include <vector>
#include <string>
#include <cstdlib>
#include <string>
#include <utility>
#include <iostream>
// root classes
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TDirectory.h"
#include "TObjString.h"
#include "TGraphErrors.h"
// tmva classes
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"

// make common namespaces implicit
using namespace std;
using namespace TMVA;

// global constants
static const UInt_t NTxt(3);
static const UInt_t NVtx(4);
static const UInt_t NHist(4);
static const UInt_t NRange(2);
//static const UInt_t NMethods(3);
static const UInt_t NMethods(1);
static const UInt_t NEneBins(10);
static const UInt_t NTmvaVar(26);
static const UInt_t NTmvaSpec(1);
static const UInt_t NCalibBins(10);

// tmva constants
static const UInt_t NTmvaHistMax(100);
static const string STmvaPrefix("TMVARegression");

// default arguments
static const string SInDef("./eicrecon_output/merged/forLowThresholdCheck.withDDSim.epic23080image.e220th45n250Kpim.d18m9y2023.plugin.root");
static const string SOutDef("forLowTresholdCheck.withDDSim.epic23080image.e220th45n250Kpim.d18m9y2023.tmva.root");
static const string STupleDef("JCalibrateHCalWithImaging/ntForCalibration");



// train and apply bhcal calibration ------------------------------------------

void TrainAndApplyBHCalCalibration(const string sInput = SInDef, const string sOutput = SOutDef, const string sTuple = STupleDef) {

  // lower verbosity
  gErrorIgnoreLevel = kWarning;
  cout << "\n  Beginning BHCal calibration training and evaluation script..." << endl;

  // options ------------------------------------------------------------------

  // general tmva parameters
  const bool   addSpectators(false);
  const float  treeWeight(1.0);
  const string sTarget("ePar");
  const string sLoader("LowThresholdCheck_DDSimPiMinus");
  const TCut   trainCut("eSumBHCal>0");

  // tmva training & spectator variables
  const string sTmvaVar[NTmvaVar] = {
    "eLeadBHCal",
    "eLeadBEMC",
    "hLeadBHCal",
    "hLeadBEMC",
    "fLeadBHCal",
    "fLeadBEMC",
    "nHitsLeadBHCal",
    "nHitsLeadBEMC",
    "eSumImage",
    "eSumSciFi",
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
    //"eSumImageLayer2",
    "eSumImageLayer3",
    "eSumImageLayer4",
    //"eSumImageLayer5",
    "eSumImageLayer6"
  };
  const string sTmvaSpec[NTmvaSpec] = {""};
  //const string sMethods[NMethods]   = {"LD", "MLP", "BDTG"};
  const string sMethods[NMethods]   = {"LD"};

  // ecal cut parameters
  const bool   doECalCut(false);
  const double eneECalRange[NRange] = {0.5, 100.};

  // histogram parameters
  const bool   isCalibrated[NHist] = {false, false, true, true};
  const string sHCalEne[NEneBins]  = {
    "hHCalEne_ene2",
    "hHCalEne_ene3",
    "hHCalEne_ene4",
    "hHCalEne_ene5",
    "hHCalEne_ene6",
    "hHCalEne_ene8",
    "hHCalEne_ene10",
    "hHCalEne_ene12",
    "hHCalEne_ene16",
    "hHCalEne_ene20"
  };
  const string sHCalDiff[NEneBins] = {
    "hHCalDiff_ene2",
    "hHCalDiff_ene3",
    "hHCalDiff_ene4",
    "hHCalDiff_ene5",
    "hHCalDiff_ene6",
    "hHCalDiff_ene8",
    "hHCalDiff_ene10",
    "hHCalDiff_ene12",
    "hHCalDiff_ene16",
    "hHCalDiff_ene20"
  };
  const string sHCalCalibBase[NCalibBins] = {
    "hHCalCalib_ene2",
    "hHCalCalib_ene3",
    "hHCalCalib_ene4",
    "hHCalCalib_ene5",
    "hHCalCalib_ene6",
    "hHCalCalib_ene8",
    "hHCalCalib_ene10",
    "hHCalCalib_ene12",
    "hHCalCalib_ene16",
    "hHCalCalib_ene20"
  };

  // generic resolution parameters
  const double enePar[NEneBins]        = {2.,  3.,  4.,  5.,  6.,  8,   10.,  12.,  16.,  20.};
  const double eneParMin[NEneBins]     = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 9.5,  11.5, 13.5, 18.5};
  const double eneParMax[NEneBins]     = {2.5, 3.5, 4.5, 5.5, 6.5, 9.5, 11.5, 13.5, 18.5, 21.5};
  const double eneCalibMin[NCalibBins] = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 9.5,  11.5, 13.5, 18.5};
  const double eneCalibMax[NCalibBins] = {2.5, 3.5, 4.5, 5.5, 6.5, 9.5, 11.5, 13.5, 18.5, 21.5};

  // reco vs. par ene resolution parameters
  const double xFitEneMin[NEneBins]  = {0., 0., 0., 1., 1.,  2.,  2.,  4.,  4.,  8.};
  const double xFitEneMax[NEneBins]  = {4., 6., 8., 9., 11., 14., 18., 20., 28., 32.};
  const double ampEneGuess[NEneBins] = {1., 1., 1., 1., 1.,  1.,  1.,  1.,  1.,  1.};
  const double muEneGuess[NEneBins]  = {2., 3., 4., 5., 6.,  8.,  10., 12., 16., 20.};
  const double sigEneGuess[NEneBins] = {1., 1., 1., 1., 1.,  1.,  3.,  3.,  3.,  7.};
  const string sFitEne[NEneBins]     = {
    "fFitEne_ene2",
    "fFitEne_ene3",
    "fFitEne_ene4",
    "fFitEne_ene5",
    "fFitEne_ene6",
    "fFitEne_ene8",
    "fFitEne_ene10",
    "fFitEne_ene12",
    "fFitEne_ene16",
    "fFitEne_ene20"
  };

  // diff vs. par ene resolution parameters
  const double xFitDiffMin[NEneBins]  = {-1., -1., -1., -1., -1., -1., -1., -1., -1., -1.};
  const double xFitDiffMax[NEneBins]  = {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.};
  const double ampDiffGuess[NEneBins] = {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.};
  const double muDiffGuess[NEneBins]  = {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.};
  const double sigDiffGuess[NEneBins] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  const string sFitDiff[NEneBins]     = {
    "fFitDiff_ene2",
    "fFitDiff_ene3",
    "fFitDiff_ene4",
    "fFitDiff_ene5",
    "fFitDiff_ene6",
    "fFitDiff_ene8",
    "fFitDiff_ene10",
    "fFitDiff_ene12",
    "fFitDiff_ene16",
    "fFitDiff_ene20"
  };

  // calibrated reco vs. par ene resolution parameters
  const double xFitCalibMin[NCalibBins]  = {0., 0., 0., 1., 1.,  2.,  2.,  4.,  4.,  8.};
  const double xFitCalibMax[NCalibBins]  = {4., 6., 8., 9., 11., 14., 18., 20., 28., 32.};
  const double ampCalibGuess[NCalibBins] = {1., 1., 1., 1., 1.,  1.,  1.,  1.,  1.,  1.};
  const double muCalibGuess[NCalibBins]  = {2., 3., 4., 5., 6.,  8.,  10., 12., 16., 20.};
  const double sigCalibGuess[NCalibBins] = {1., 1., 1., 1., 1.,  1.,  3.,  3.,  3.,  7.};
  const string sFitCalibBase[NCalibBins] = {
    "fFitCalib_ene2",
    "fFitCalib_ene3",
    "fFitCalib_ene4",
    "fFitCalib_ene5",
    "fFitCalib_ene6",
    "fFitCalib_ene8",
    "fFitCalib_ene10",
    "fFitCalib_ene12",
    "fFitCalib_ene16",
    "fFitCalib_ene20"
  };

  // load input ---------------------------------------------------------------

  // open files
  TFile *fOutput  = new TFile(sOutput.data(), "recreate");
  TFile *fInTrain = new TFile(sInput.data(),  "read");
  TFile *fInApply = new TFile(sInput.data(),  "read");
  if (!fOutput || !fInTrain || !fInApply) {
    cerr << "PANIC: couldn't open a file!\n"
         << "        fInTrain = " << fInTrain << ", fInApply = " << fInApply << ", fOutput = " << fOutput << "\n"
         << endl;
    return;
  }
  cout << "    Opened files:\n"
       << "      fInput  = " << sInput.data() << "\n"
       << "      fOutput = " << sOutput.data()
       << endl;

  // grab input tuple
  TNtuple *ntToTrain = (TNtuple*) fInTrain -> Get(sTuple.data());
  TNtuple *ntToApply = (TNtuple*) fInApply -> Get(sTuple.data());
  if (!ntToTrain || !ntToApply) {
    cerr << "PANIC: couldn't grab input tuple!\n"
         << "       name  = " << sTuple    << "\n"
         << "       train = " << ntToTrain << "\n"
         << "       apply = " << ntToApply << "\n"
         << endl;
    return;
  }
  cout << "    Grabbed input tuples:\n"
       << "      tuple = " << sTuple
       << endl;

  // declare tuple leaves
  float ePar_train;
  float fracParVsLeadBHCal_train;
  float fracParVsLeadBEMC_train;
  float fracParVsSumBHCal_train;
  float fracParVsSumBEMC_train;
  float fracLeadBHCalVsBEMC_train;
  float fracSumBHCalVsBEMC_train;
  float eLeadBHCal_train;
  float eLeadBEMC_train;
  float eSumBHCal_train;
  float eSumBEMC_train;
  float diffLeadBHCal_train;
  float diffLeadBEMC_train;
  float diffSumBHCal_train;
  float diffSumBEMC_train;
  float nHitsLeadBHCal_train;
  float nHitsLeadBEMC_train;
  float nClustBHCal_train;
  float nClustBEMC_train;
  float hLeadBHCal_train;
  float hLeadBEMC_train;
  float fLeadBHCal_train;
  float fLeadBEMC_train;
  float eLeadImage_train;
  float eSumImage_train;
  float eLeadSciFi_train;
  float eSumSciFi_train;
  float nClustImage_train;
  float nClustSciFi_train;
  float hLeadImage_train;
  float hLeadSciFi_train;
  float fLeadImage_train;
  float fLeadSciFi_train;
  float eSumSciFiLayer1_train;
  float eSumSciFiLayer2_train;
  float eSumSciFiLayer3_train;
  float eSumSciFiLayer4_train;
  float eSumSciFiLayer5_train;
  float eSumSciFiLayer6_train;
  float eSumSciFiLayer7_train;
  float eSumSciFiLayer8_train;
  float eSumSciFiLayer9_train;
  float eSumSciFiLayer10_train;
  float eSumSciFiLayer11_train;
  float eSumSciFiLayer12_train;
  float eSumImageLayer1_train;
  float eSumImageLayer2_train;
  float eSumImageLayer3_train;
  float eSumImageLayer4_train;
  float eSumImageLayer5_train;
  float eSumImageLayer6_train;

  float ePar_apply;
  float fracParVsLeadBHCal_apply;
  float fracParVsLeadBEMC_apply;
  float fracParVsSumBHCal_apply;
  float fracParVsSumBEMC_apply;
  float fracLeadBHCalVsBEMC_apply;
  float fracSumBHCalVsBEMC_apply;
  float eLeadBHCal_apply;
  float eLeadBEMC_apply;
  float eSumBHCal_apply;
  float eSumBEMC_apply;
  float diffLeadBHCal_apply;
  float diffLeadBEMC_apply;
  float diffSumBHCal_apply;
  float diffSumBEMC_apply;
  float nHitsLeadBHCal_apply;
  float nHitsLeadBEMC_apply;
  float nClustBHCal_apply;
  float nClustBEMC_apply;
  float hLeadBHCal_apply;
  float hLeadBEMC_apply;
  float fLeadBHCal_apply;
  float fLeadBEMC_apply;
  float eLeadImage_apply;
  float eSumImage_apply;
  float eLeadSciFi_apply;
  float eSumSciFi_apply;
  float nClustImage_apply;
  float nClustSciFi_apply;
  float hLeadImage_apply;
  float hLeadSciFi_apply;
  float fLeadImage_apply;
  float fLeadSciFi_apply;
  float eSumSciFiLayer1_apply;
  float eSumSciFiLayer2_apply;
  float eSumSciFiLayer3_apply;
  float eSumSciFiLayer4_apply;
  float eSumSciFiLayer5_apply;
  float eSumSciFiLayer6_apply;
  float eSumSciFiLayer7_apply;
  float eSumSciFiLayer8_apply;
  float eSumSciFiLayer9_apply;
  float eSumSciFiLayer10_apply;
  float eSumSciFiLayer11_apply;
  float eSumSciFiLayer12_apply;
  float eSumImageLayer1_apply;
  float eSumImageLayer2_apply;
  float eSumImageLayer3_apply;
  float eSumImageLayer4_apply;
  float eSumImageLayer5_apply;
  float eSumImageLayer6_apply;

  // set tuple branches
  ntToTrain -> SetBranchAddress("ePar",                &ePar_train);
  ntToTrain -> SetBranchAddress("fracParVsLeadBHCal",  &fracParVsLeadBHCal_train);
  ntToTrain -> SetBranchAddress("fracParVsLeadBEMC",   &fracParVsLeadBEMC_train);
  ntToTrain -> SetBranchAddress("fracParVsSumBHCal",   &fracParVsSumBHCal_train);
  ntToTrain -> SetBranchAddress("fracParVsSumBEMC",    &fracParVsSumBEMC_train);
  ntToTrain -> SetBranchAddress("fracLeadBHCalVsBEMC", &fracLeadBHCalVsBEMC_train);
  ntToTrain -> SetBranchAddress("fracSumBHCalVsBEMC",  &fracSumBHCalVsBEMC_train);
  ntToTrain -> SetBranchAddress("eLeadBHCal",          &eLeadBHCal_train);
  ntToTrain -> SetBranchAddress("eLeadBEMC",           &eLeadBEMC_train);
  ntToTrain -> SetBranchAddress("eSumBHCal",           &eSumBHCal_train);
  ntToTrain -> SetBranchAddress("eSumBEMC",            &eSumBEMC_train);
  ntToTrain -> SetBranchAddress("diffLeadBHCal",       &diffLeadBHCal_train);
  ntToTrain -> SetBranchAddress("diffLeadBEMC",        &diffLeadBEMC_train);
  ntToTrain -> SetBranchAddress("diffSumBHCal",        &diffSumBHCal_train);
  ntToTrain -> SetBranchAddress("diffSumBEMC",         &diffSumBEMC_train);
  ntToTrain -> SetBranchAddress("nHitsLeadBHCal",      &nHitsLeadBHCal_train);
  ntToTrain -> SetBranchAddress("nHitsLeadBEMC",       &nHitsLeadBEMC_train);
  ntToTrain -> SetBranchAddress("nClustBHCal",         &nClustBHCal_train);
  ntToTrain -> SetBranchAddress("nClustBEMC",          &nClustBEMC_train);
  ntToTrain -> SetBranchAddress("hLeadBHCal",          &hLeadBHCal_train);
  ntToTrain -> SetBranchAddress("hLeadBEMC",           &hLeadBEMC_train);
  ntToTrain -> SetBranchAddress("fLeadBHCal",          &fLeadBHCal_train);
  ntToTrain -> SetBranchAddress("fLeadBEMC",           &fLeadBEMC_train);
  ntToTrain -> SetBranchAddress("eLeadImage",          &eLeadImage_train);
  ntToTrain -> SetBranchAddress("eSumImage",           &eSumImage_train);
  ntToTrain -> SetBranchAddress("eLeadSciFi",          &eLeadSciFi_train);
  ntToTrain -> SetBranchAddress("eSumSciFi",           &eSumSciFi_train);
  ntToTrain -> SetBranchAddress("nClustImage",         &nClustImage_train);
  ntToTrain -> SetBranchAddress("nClustSciFi",         &nClustSciFi_train);
  ntToTrain -> SetBranchAddress("hLeadImage",          &hLeadImage_train);
  ntToTrain -> SetBranchAddress("hLeadSciFi",          &hLeadSciFi_train);
  ntToTrain -> SetBranchAddress("fLeadImage",          &fLeadImage_train);
  ntToTrain -> SetBranchAddress("fLeadSciFi",          &fLeadSciFi_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer1",     &eSumSciFiLayer1_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer2",     &eSumSciFiLayer2_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer3",     &eSumSciFiLayer3_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer4",     &eSumSciFiLayer4_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer5",     &eSumSciFiLayer5_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer6",     &eSumSciFiLayer6_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer7",     &eSumSciFiLayer7_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer8",     &eSumSciFiLayer8_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer9",     &eSumSciFiLayer9_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer10",    &eSumSciFiLayer10_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer11",    &eSumSciFiLayer11_train);
  ntToTrain -> SetBranchAddress("eSumSciFiLayer12",    &eSumSciFiLayer12_train);
  ntToTrain -> SetBranchAddress("eSumImageLayer1",     &eSumImageLayer1_train);
  ntToTrain -> SetBranchAddress("eSumImageLayer2",     &eSumImageLayer2_train);
  ntToTrain -> SetBranchAddress("eSumImageLayer3",     &eSumImageLayer3_train);
  ntToTrain -> SetBranchAddress("eSumImageLayer4",     &eSumImageLayer4_train);
  ntToTrain -> SetBranchAddress("eSumImageLayer5",     &eSumImageLayer5_train);
  ntToTrain -> SetBranchAddress("eSumImageLayer6",     &eSumImageLayer6_train);

  ntToApply -> SetBranchAddress("ePar",                &ePar_apply);
  ntToApply -> SetBranchAddress("fracParVsLeadBHCal",  &fracParVsLeadBHCal_apply);
  ntToApply -> SetBranchAddress("fracParVsLeadBEMC",   &fracParVsLeadBEMC_apply);
  ntToApply -> SetBranchAddress("fracParVsSumBHCal",   &fracParVsSumBHCal_apply);
  ntToApply -> SetBranchAddress("fracParVsSumBEMC",    &fracParVsSumBEMC_apply);
  ntToApply -> SetBranchAddress("fracLeadBHCalVsBEMC", &fracLeadBHCalVsBEMC_apply);
  ntToApply -> SetBranchAddress("fracSumBHCalVsBEMC",  &fracSumBHCalVsBEMC_apply);
  ntToApply -> SetBranchAddress("eLeadBHCal",          &eLeadBHCal_apply);
  ntToApply -> SetBranchAddress("eLeadBEMC",           &eLeadBEMC_apply);
  ntToApply -> SetBranchAddress("eSumBHCal",           &eSumBHCal_apply);
  ntToApply -> SetBranchAddress("eSumBEMC",            &eSumBEMC_apply);
  ntToApply -> SetBranchAddress("diffLeadBHCal",       &diffLeadBHCal_apply);
  ntToApply -> SetBranchAddress("diffLeadBEMC",        &diffLeadBEMC_apply);
  ntToApply -> SetBranchAddress("diffSumBHCal",        &diffSumBHCal_apply);
  ntToApply -> SetBranchAddress("diffSumBEMC",         &diffSumBEMC_apply);
  ntToApply -> SetBranchAddress("nHitsLeadBHCal",      &nHitsLeadBHCal_apply);
  ntToApply -> SetBranchAddress("nHitsLeadBEMC",       &nHitsLeadBEMC_apply);
  ntToApply -> SetBranchAddress("nClustBHCal",         &nClustBHCal_apply);
  ntToApply -> SetBranchAddress("nClustBEMC",          &nClustBEMC_apply);
  ntToApply -> SetBranchAddress("hLeadBHCal",          &hLeadBHCal_apply);
  ntToApply -> SetBranchAddress("hLeadBEMC",           &hLeadBEMC_apply);
  ntToApply -> SetBranchAddress("fLeadBHCal",          &fLeadBHCal_apply);
  ntToApply -> SetBranchAddress("fLeadBEMC",           &fLeadBEMC_apply);
  ntToApply -> SetBranchAddress("eLeadImage",          &eLeadImage_apply);
  ntToApply -> SetBranchAddress("eSumImage",           &eSumImage_apply);
  ntToApply -> SetBranchAddress("eLeadSciFi",          &eLeadSciFi_apply);
  ntToApply -> SetBranchAddress("eSumSciFi",           &eSumSciFi_apply);
  ntToApply -> SetBranchAddress("nClustImage",         &nClustImage_apply);
  ntToApply -> SetBranchAddress("nClustSciFi",         &nClustSciFi_apply);
  ntToApply -> SetBranchAddress("hLeadImage",          &hLeadImage_apply);
  ntToApply -> SetBranchAddress("hLeadSciFi",          &hLeadSciFi_apply);
  ntToApply -> SetBranchAddress("fLeadImage",          &fLeadImage_apply);
  ntToApply -> SetBranchAddress("fLeadSciFi",          &fLeadSciFi_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer1",     &eSumSciFiLayer1_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer2",     &eSumSciFiLayer2_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer3",     &eSumSciFiLayer3_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer4",     &eSumSciFiLayer4_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer5",     &eSumSciFiLayer5_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer6",     &eSumSciFiLayer6_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer7",     &eSumSciFiLayer7_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer8",     &eSumSciFiLayer8_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer9",     &eSumSciFiLayer9_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer10",    &eSumSciFiLayer10_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer11",    &eSumSciFiLayer11_apply);
  ntToApply -> SetBranchAddress("eSumSciFiLayer12",    &eSumSciFiLayer12_apply);
  ntToApply -> SetBranchAddress("eSumImageLayer1",     &eSumImageLayer1_apply);
  ntToApply -> SetBranchAddress("eSumImageLayer2",     &eSumImageLayer2_apply);
  ntToApply -> SetBranchAddress("eSumImageLayer3",     &eSumImageLayer3_apply);
  ntToApply -> SetBranchAddress("eSumImageLayer4",     &eSumImageLayer4_apply);
  ntToApply -> SetBranchAddress("eSumImageLayer5",     &eSumImageLayer5_apply);
  ntToApply -> SetBranchAddress("eSumImageLayer6",     &eSumImageLayer6_apply);
  cout << "    Set tuple branches." << endl;

  // declare output histograms ------------------------------------------------

  // general histograms & profiles
  TH1D     *hHCalFrac[NHist];
  TH1D     *hHCalDiff[NHist];
  TH1D     *hECalFrac[NHist];
  TH1D     *hECalDiff[NHist];
  TH2D     *hHCalEneVsPar[NHist];
  TH2D     *hECalEneVsPar[NHist];
  TH2D     *hHCalFracVsPar[NHist];
  TH2D     *hHCalDiffVsPar[NHist];
  TH2D     *hECalFracVsPar[NHist];
  TH2D     *hECalDiffVsPar[NHist];
  TH2D     *hHCalVsECalFrac[NHist];
  TH2D     *hHCalVsECalDiff[NHist];
  TH2D     *hHCalFracVsTotalFrac[NHist];
  TH2D     *hHCalDiffVsTotalFrac[NHist];
  TH2D     *hECalFracVsTotalFrac[NHist];
  TH2D     *hECalDiffVsTotalFrac[NHist];
  TProfile *pHCalEneVsPar[NHist];
  TProfile *pECalEneVsPar[NHist];
  TProfile *pHCalFracVsPar[NHist];
  TProfile *pHCalDiffVsPar[NHist];
  TProfile *pECalFracVsPar[NHist];
  TProfile *pECalDiffVsPar[NHist];
  TProfile *pHCalVsECalFrac[NHist];
  TProfile *pHCalVsECalDiff[NHist];
  TProfile *pHCalFracVsTotalFrac[NHist];
  TProfile *pHCalDiffVsTotalFrac[NHist];
  TProfile *pECalFracVsTotalFrac[NHist];
  TProfile *pECalDiffVsTotalFrac[NHist];

  // resolution histograms
  TH1D *hHCalEneBin[NEneBins];
  TH1D *hHCalDiffBin[NEneBins];

  // histogram binning
  const size_t nEneBins(41);
  const size_t nEneBins2D(410);
  const size_t nDiffBins(700);
  const size_t nFracBins(305);
  const float  rEneBins[NRange]  = {-1.,   40.};
  const float  rDiffBins[NRange] = {-1.5,  5.5};
  const float  rFracBins[NRange] = {-0.05, 3.};

  // declare uncalibrated histograms
  hHCalFrac[0]            = new TH1D("hLeadHCalFrac_uncal",            "", nFracBins, rFracBins[0], rFracBins[1]);
  hHCalFrac[1]            = new TH1D("hSumHCalFrac_uncal",             "", nFracBins, rFracBins[0], rFracBins[1]);
  hHCalDiff[0]            = new TH1D("hLeadHCalDiff_uncal",            "", nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalDiff[1]            = new TH1D("hSumHCalDiff_uncal",             "", nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalFrac[0]            = new TH1D("hLeadECalFrac_uncal",            "", nFracBins, rFracBins[0], rFracBins[1]);
  hECalFrac[1]            = new TH1D("hSumECalFrac_uncal",             "", nFracBins, rFracBins[0], rFracBins[1]);
  hECalDiff[0]            = new TH1D("hLeadECalDiff_uncal",            "", nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalDiff[1]            = new TH1D("hSumECalDiff_uncal",             "", nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalEneVsPar[0]        = new TH2D("hLeadHCalVsParEne_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nEneBins,  rEneBins[0],  rEneBins[1]);
  hHCalEneVsPar[1]        = new TH2D("hSumHCalVsParEne_uncal",         "", nEneBins,  rEneBins[0],  rEneBins[1],  nEneBins,  rEneBins[0],  rEneBins[1]);     
  hECalEneVsPar[0]        = new TH2D("hLeadECalVsParEne_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nEneBins,  rEneBins[0],  rEneBins[1]);
  hECalEneVsPar[1]        = new TH2D("hSumECalVsParEne_uncal",         "", nEneBins,  rEneBins[0],  rEneBins[1],  nEneBins,  rEneBins[0],  rEneBins[1]);     
  hHCalFracVsPar[0]       = new TH2D("hLeadHCalFracVsPar_uncal",       "", nEneBins,  rEneBins[0],  rEneBins[1],  nFracBins, rFracBins[0], rFracBins[1]);
  hHCalFracVsPar[1]       = new TH2D("hSumHCalFracVsPar_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nFracBins, rFracBins[0], rFracBins[1]);
  hHCalDiffVsPar[0]       = new TH2D("hLeadHCalDiffVsPar_uncal",       "", nEneBins,  rEneBins[0],  rEneBins[1],  nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalDiffVsPar[1]       = new TH2D("hSumHCalDiffVsPar_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalFracVsPar[0]       = new TH2D("hLeadECalFracVsPar_uncal",       "", nEneBins,  rEneBins[0],  rEneBins[1],  nFracBins, rFracBins[0], rFracBins[1]);
  hECalFracVsPar[1]       = new TH2D("hSumECalFracVsPar_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nFracBins, rFracBins[0], rFracBins[1]);
  hECalDiffVsPar[0]       = new TH2D("hLeadECalDiffVsPar_uncal",       "", nEneBins,  rEneBins[0],  rEneBins[1],  nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalDiffVsPar[1]       = new TH2D("hSumECalDiffVsPar_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalVsECalFrac[0]      = new TH2D("hLeadHCalVsLeadECalFrac_uncal",  "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hHCalVsECalFrac[1]      = new TH2D("hSumHCalVsSumECalFrac_uncal",    "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hHCalVsECalDiff[0]      = new TH2D("hLeadHCalVsLeadECalDiff_uncal",  "", nDiffBins, rDiffBins[0], rDiffBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalVsECalDiff[1]      = new TH2D("hSumHCalVsSumECalDiff_uncal",    "", nDiffBins, rDiffBins[0], rDiffBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalFracVsTotalFrac[0] = new TH2D("hLeadHCalFracVsTotalFrac_uncal", "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hHCalFracVsTotalFrac[1] = new TH2D("hSumHCalFracVsTotalFrac_uncal",  "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hHCalDiffVsTotalFrac[0] = new TH2D("hLeadHCalDiffVsTotalFrac_uncal", "", nFracBins, rFracBins[0], rFracBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalDiffVsTotalFrac[1] = new TH2D("hSumHCalDiffVsTotalFrac_uncal",  "", nFracBins, rFracBins[0], rFracBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalFracVsTotalFrac[0] = new TH2D("hLeadECalFracVsTotalFrac_uncal", "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hECalFracVsTotalFrac[1] = new TH2D("hSumECalFracVsTotalFrac_uncal",  "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hECalDiffVsTotalFrac[0] = new TH2D("hLeadECalDiffVsTotalFrac_uncal", "", nFracBins, rFracBins[0], rFracBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalDiffVsTotalFrac[1] = new TH2D("hSumECalDiffVsTotalFrac_uncal",  "", nFracBins, rFracBins[0], rFracBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  // declare calibrated histograms
  hHCalFrac[2]            = new TH1D("hLeadHCalFrac_calib",            "", nFracBins, rFracBins[0], rFracBins[1]);
  hHCalFrac[3]            = new TH1D("hSumHCalFrac_calib",             "", nFracBins, rFracBins[0], rFracBins[1]);
  hHCalDiff[2]            = new TH1D("hLeadHCalDiff_calib",            "", nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalDiff[3]            = new TH1D("hSumHCalDiff_calib",             "", nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalFrac[2]            = new TH1D("hLeadECalFrac_calib",            "", nFracBins, rFracBins[0], rFracBins[1]);
  hECalFrac[3]            = new TH1D("hSumECalFrac_calib",             "", nFracBins, rFracBins[0], rFracBins[1]);
  hECalDiff[2]            = new TH1D("hLeadECalDiff_calib",            "", nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalDiff[3]            = new TH1D("hSumECalDiff_calib",             "", nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalEneVsPar[2]        = new TH2D("hLeadHCalVsParEne_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nEneBins,  rEneBins[0],  rEneBins[1]);
  hHCalEneVsPar[3]        = new TH2D("hSumHCalVsParEne_calib",         "", nEneBins,  rEneBins[0],  rEneBins[1],  nEneBins,  rEneBins[0],  rEneBins[1]);     
  hECalEneVsPar[2]        = new TH2D("hLeadECalVsParEne_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nEneBins,  rEneBins[0],  rEneBins[1]);
  hECalEneVsPar[3]        = new TH2D("hSumECalVsParEne_calib",         "", nEneBins,  rEneBins[0],  rEneBins[1],  nEneBins,  rEneBins[0],  rEneBins[1]);     
  hHCalFracVsPar[2]       = new TH2D("hLeadHCalFracVsPar_calib",       "", nEneBins,  rEneBins[0],  rEneBins[1],  nFracBins, rFracBins[0], rFracBins[1]);
  hHCalFracVsPar[3]       = new TH2D("hSumHCalFracVsPar_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nFracBins, rFracBins[0], rFracBins[1]);
  hHCalDiffVsPar[2]       = new TH2D("hLeadHCalDiffVsPar_calib",       "", nEneBins,  rEneBins[0],  rEneBins[1],  nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalDiffVsPar[3]       = new TH2D("hSumHCalDiffVsPar_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalFracVsPar[2]       = new TH2D("hLeadECalFracVsPar_calib",       "", nEneBins,  rEneBins[0],  rEneBins[1],  nFracBins, rFracBins[0], rFracBins[1]);
  hECalFracVsPar[3]       = new TH2D("hSumECalFracVsPar_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nFracBins, rFracBins[0], rFracBins[1]);
  hECalDiffVsPar[2]       = new TH2D("hLeadECalDiffVsPar_calib",       "", nEneBins,  rEneBins[0],  rEneBins[1],  nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalDiffVsPar[3]       = new TH2D("hSumECalDiffVsPar_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalVsECalFrac[2]      = new TH2D("hLeadHCalVsLeadECalFrac_calib",  "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hHCalVsECalFrac[3]      = new TH2D("hSumHCalVsSumECalFrac_calib",    "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hHCalVsECalDiff[2]      = new TH2D("hLeadHCalVsLeadECalDiff_calib",  "", nDiffBins, rDiffBins[0], rDiffBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalVsECalDiff[3]      = new TH2D("hSumHCalVsSumECalDiff_calib",    "", nDiffBins, rDiffBins[0], rDiffBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalFracVsTotalFrac[2] = new TH2D("hLeadHCalFracVsTotalFrac_calib", "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hHCalFracVsTotalFrac[3] = new TH2D("hSumHCalFracVsTotalFrac_calib",  "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hHCalDiffVsTotalFrac[2] = new TH2D("hLeadHCalDiffVsTotalFrac_calib", "", nFracBins, rFracBins[0], rFracBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  hHCalDiffVsTotalFrac[3] = new TH2D("hSumHCalDiffVsTotalFrac_calib",  "", nFracBins, rFracBins[0], rFracBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalFracVsTotalFrac[2] = new TH2D("hLeadECalFracVsTotalFrac_calib", "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hECalFracVsTotalFrac[3] = new TH2D("hSumECalFracVsTotalFrac_calib",  "", nFracBins, rFracBins[0], rFracBins[1], nFracBins, rFracBins[0], rFracBins[1]);
  hECalDiffVsTotalFrac[2] = new TH2D("hLeadECalDiffVsTotalFrac_calib", "", nFracBins, rFracBins[0], rFracBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  hECalDiffVsTotalFrac[3] = new TH2D("hSumECalDiffVsTotalFrac_calib",  "", nFracBins, rFracBins[0], rFracBins[1], nDiffBins, rDiffBins[0], rDiffBins[1]);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hHCalFrac[iHist]            -> Sumw2();
    hHCalDiff[iHist]            -> Sumw2();
    hECalFrac[iHist]            -> Sumw2();
    hECalDiff[iHist]            -> Sumw2();
    hHCalEneVsPar[iHist]        -> Sumw2();
    hECalEneVsPar[iHist]        -> Sumw2();
    hHCalFracVsPar[iHist]       -> Sumw2();
    hHCalDiffVsPar[iHist]       -> Sumw2();
    hECalFracVsPar[iHist]       -> Sumw2();
    hECalDiffVsPar[iHist]       -> Sumw2();
    hHCalVsECalFrac[iHist]      -> Sumw2();
    hHCalVsECalDiff[iHist]      -> Sumw2();
    hHCalFracVsTotalFrac[iHist] -> Sumw2();
    hHCalDiffVsTotalFrac[iHist] -> Sumw2();
    hECalFracVsTotalFrac[iHist] -> Sumw2();
    hECalDiffVsTotalFrac[iHist] -> Sumw2();
  }

  // declare uncalibrated profiles
  pHCalEneVsPar[0]        = new TProfile("pLeadHCalVsParEne_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pHCalEneVsPar[1]        = new TProfile("pSumHCalVsParEne_uncal",         "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");     
  pECalEneVsPar[0]        = new TProfile("pLeadECalVsParEne_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pECalEneVsPar[1]        = new TProfile("pSumECalVsParEne_uncal",         "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");     
  pHCalFracVsPar[0]       = new TProfile("pLeadHCalFracVsPar_uncal",       "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pHCalFracVsPar[1]       = new TProfile("pSumHCalFracVsPar_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pHCalDiffVsPar[0]       = new TProfile("pLeadHCalDiffVsPar_uncal",       "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pHCalDiffVsPar[1]       = new TProfile("pSumHCalDiffVsPar_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pECalFracVsPar[0]       = new TProfile("pLeadECalFracVsPar_uncal",       "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pECalFracVsPar[1]       = new TProfile("pSumECalFracVsPar_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pECalDiffVsPar[0]       = new TProfile("pLeadECalDiffVsPar_uncal",       "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pECalDiffVsPar[1]       = new TProfile("pSumECalDiffVsPar_uncal",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pHCalVsECalFrac[0]      = new TProfile("pLeadHCalVsLeadECalFrac_uncal",  "", nFracBins, rFracBins[0], rFracBins[1], "S");
  pHCalVsECalFrac[1]      = new TProfile("pSumHCalVsSumECalFrac_uncal",    "", nFracBins, rFracBins[0], rFracBins[1], "S");
  pHCalVsECalDiff[0]      = new TProfile("pLeadHCalVsLeadECalDiff_uncal",  "", nDiffBins, rDiffBins[0], rDiffBins[1], "S");
  pHCalVsECalDiff[1]      = new TProfile("pSumHCalVsSumECalDiff_uncal",    "", nDiffBins, rDiffBins[0], rDiffBins[1], "S");
  pHCalFracVsTotalFrac[0] = new TProfile("pLeadHCalFracVsTotalFrac_uncal", "", nFracBins, rFracBins[0], rFracBins[1], "S"); 
  pHCalFracVsTotalFrac[1] = new TProfile("pSumHCalFracVsTotalFrac_uncal",  "", nFracBins, rFracBins[0], rFracBins[1], "S");
  pHCalDiffVsTotalFrac[0] = new TProfile("pLeadHCalDiffVsTotalFrac_uncal", "", nFracBins, rFracBins[0], rFracBins[1], "S"); 
  pHCalDiffVsTotalFrac[1] = new TProfile("pSumHCalDiffVsTotalFrac_uncal",  "", nDiffBins, rDiffBins[0], rDiffBins[1], "S");
  pECalFracVsTotalFrac[0] = new TProfile("pLeadECalFracVsTotalFrac_uncal", "", nFracBins, rFracBins[0], rFracBins[1], "S"); 
  pECalFracVsTotalFrac[1] = new TProfile("pSumECalFracVsTotalFrac_uncal",  "", nFracBins, rFracBins[0], rFracBins[1], "S");
  pECalDiffVsTotalFrac[0] = new TProfile("pLeadECalDiffVsTotalFrac_uncal", "", nFracBins, rFracBins[0], rFracBins[1], "S"); 
  pECalDiffVsTotalFrac[1] = new TProfile("pSumECalDiffVsTotalFrac_uncal",  "", nDiffBins, rDiffBins[0], rDiffBins[1], "S");
  // declare calibrated profiles
  pHCalEneVsPar[2]        = new TProfile("pLeadHCalVsParEne_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pHCalEneVsPar[3]        = new TProfile("pSumHCalVsParEne_calib",         "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");     
  pECalEneVsPar[2]        = new TProfile("pLeadECalVsParEne_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pECalEneVsPar[3]        = new TProfile("pSumECalVsParEne_calib",         "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");     
  pHCalFracVsPar[2]       = new TProfile("pLeadHCalFracVsPar_calib",       "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pHCalFracVsPar[3]       = new TProfile("pSumHCalFracVsPar_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pHCalDiffVsPar[2]       = new TProfile("pLeadHCalDiffVsPar_calib",       "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pHCalDiffVsPar[3]       = new TProfile("pSumHCalDiffVsPar_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pECalFracVsPar[2]       = new TProfile("pLeadECalFracVsPar_calib",       "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pECalFracVsPar[3]       = new TProfile("pSumECalFracVsPar_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pECalDiffVsPar[2]       = new TProfile("pLeadECalDiffVsPar_calib",       "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pECalDiffVsPar[3]       = new TProfile("pSumECalDiffVsPar_calib",        "", nEneBins,  rEneBins[0],  rEneBins[1],  "S");
  pHCalVsECalFrac[2]      = new TProfile("pLeadHCalVsLeadECalFrac_calib",  "", nFracBins, rFracBins[0], rFracBins[1], "S");
  pHCalVsECalFrac[3]      = new TProfile("pSumHCalVsSumECalFrac_calib",    "", nFracBins, rFracBins[0], rFracBins[1], "S");
  pHCalVsECalDiff[2]      = new TProfile("pLeadHCalVsLeadECalDiff_calib",  "", nDiffBins, rDiffBins[0], rDiffBins[1], "S");
  pHCalVsECalDiff[3]      = new TProfile("pSumHCalVsSumECalDiff_calib",    "", nDiffBins, rDiffBins[0], rDiffBins[1], "S");
  pHCalFracVsTotalFrac[2] = new TProfile("pLeadHCalFracVsTotalFrac_calib", "", nFracBins, rFracBins[0], rFracBins[1], "S"); 
  pHCalFracVsTotalFrac[3] = new TProfile("pSumHCalFracVsTotalFrac_calib",  "", nFracBins, rFracBins[0], rFracBins[1], "S");
  pHCalDiffVsTotalFrac[2] = new TProfile("pLeadHCalDiffVsTotalFrac_calib", "", nFracBins, rFracBins[0], rFracBins[1], "S"); 
  pHCalDiffVsTotalFrac[3] = new TProfile("pSumHCalDiffVsTotalFrac_calib",  "", nDiffBins, rDiffBins[0], rDiffBins[1], "S");
  pECalFracVsTotalFrac[2] = new TProfile("pLeadECalFracVsTotalFrac_calib", "", nFracBins, rFracBins[0], rFracBins[1], "S"); 
  pECalFracVsTotalFrac[3] = new TProfile("pSumECalFracVsTotalFrac_calib",  "", nFracBins, rFracBins[0], rFracBins[1], "S");
  pECalDiffVsTotalFrac[2] = new TProfile("pLeadECalDiffVsTotalFrac_calib", "", nFracBins, rFracBins[0], rFracBins[1], "S"); 
  pECalDiffVsTotalFrac[3] = new TProfile("pSumECalDiffVsTotalFrac_calib",  "", nDiffBins, rDiffBins[0], rDiffBins[1], "S");

  // declare resolution histograms
  for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {
    hHCalEneBin[iEneBin]  = new TH1D(sHCalEne[iEneBin].data(),  "", nEneBins,  rEneBins[0],  rEneBins[1]);
    hHCalDiffBin[iEneBin] = new TH1D(sHCalDiff[iEneBin].data(), "", nDiffBins, rDiffBins[0], rDiffBins[1]);
    hHCalEneBin[iEneBin]  -> Sumw2();
    hHCalDiffBin[iEneBin] -> Sumw2();
  }

  // calibrated resolution histograms
  TH1D *hHCalCalibBin[NMethods][NCalibBins];
  TH2D *hCalibCalibVsPar[NMethods];
  TH2D *hHCalCalibVsPar[NMethods];
  TH2D *hHCalCalibVsCalib[NMethods];
  TH2D *hHCalCalibVsECal[NMethods];
  TH2D *hECalCalibVsPar[NMethods];
  TH2D *hECalCalibVsCalib[NMethods];

  // declare calibrated resolution histograms
  for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
    for (UInt_t iCalibBin = 0; iCalibBin < NCalibBins; iCalibBin++) {

      // make name
      TString sHCalCalib(sHCalCalibBase[iCalibBin].data());
      sHCalCalib.Append("_");
      sHCalCalib.Append(sMethods[iMethod].data());
     
      hHCalCalibBin[iMethod][iCalibBin] = new TH1D(sHCalCalib.Data(), "", nEneBins, rEneBins[0], rEneBins[1]);
      hHCalCalibBin[iMethod][iCalibBin] -> Sumw2();
    }

    // make name
    TString sCalibCalibVsPar("hCalibCalibVsPar");
    TString sHCalCalibVsPar("hHCalCalibVsPar");
    TString sHCalCalibVsCalib("hHCalCalibVsCalib");
    TString sHCalCalibVsECal("hHCalCalibVsECal");
    TString sECalCalibVsPar("hECalCalibVsPar");
    TString sECalCalibVsCalib("hECalCalibVsCalib");
    sCalibCalibVsPar.Append("_");
    sHCalCalibVsPar.Append("_");
    sHCalCalibVsCalib.Append("_");
    sHCalCalibVsECal.Append("_");
    sECalCalibVsPar.Append("_");
    sECalCalibVsCalib.Append("_");
    sCalibCalibVsPar.Append(sMethods[iMethod].data());
    sHCalCalibVsPar.Append(sMethods[iMethod].data());
    sHCalCalibVsCalib.Append(sMethods[iMethod].data());
    sHCalCalibVsECal.Append(sMethods[iMethod].data());
    sECalCalibVsPar.Append(sMethods[iMethod].data());
    sECalCalibVsCalib.Append(sMethods[iMethod].data());

    hCalibCalibVsPar[iMethod]  = new TH2D(sCalibCalibVsPar.Data(),  "", nEneBins,   rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hHCalCalibVsPar[iMethod]   = new TH2D(sHCalCalibVsPar.Data(),   "", nEneBins2D, rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hHCalCalibVsCalib[iMethod] = new TH2D(sHCalCalibVsCalib.Data(), "", nEneBins2D, rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hHCalCalibVsECal[iMethod]  = new TH2D(sHCalCalibVsECal.Data(),  "", nEneBins2D, rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hECalCalibVsPar[iMethod]   = new TH2D(sECalCalibVsPar.Data(),   "", nEneBins2D, rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hECalCalibVsCalib[iMethod] = new TH2D(sECalCalibVsCalib.Data(), "", nEneBins2D, rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hCalibCalibVsPar[iMethod]  -> Sumw2();
    hHCalCalibVsPar[iMethod]   -> Sumw2();
    hHCalCalibVsCalib[iMethod] -> Sumw2();
    hHCalCalibVsECal[iMethod]  -> Sumw2();
    hECalCalibVsPar[iMethod]   -> Sumw2();
    hECalCalibVsCalib[iMethod] -> Sumw2();
  }  // end method loop
  cout << "    Declared resolution histograms." << endl;

  // loop over ntuple entries -------------------------------------------------

  // prepare for uncalibrated tuple loop
  Long64_t nEvtsToTrain = ntToTrain -> GetEntries();
  cout << "    Looping over uncalibrated tuple: " << nEvtsToTrain << " events to process." << endl;

  Long64_t nBytes(0);
  for (Long64_t iEvtToTrain = 0; iEvtToTrain < nEvtsToTrain; iEvtToTrain++) {

    const Long64_t bytes = ntToTrain -> GetEntry(iEvtToTrain);
    if (bytes < 0.) {
      cerr << "WARNING something wrong with event " << iEvtToTrain << "! Aborting loop!" << endl;
      break;
    }
    nBytes += bytes;

    // announce progress
    const Long64_t iProgToTrain = iEvtToTrain + 1;
    if (iProgToTrain == nEvtsToTrain) {
      cout << "      Proceesing event " << iProgToTrain << "/" << nEvtsToTrain << "..." << endl;
    } else {
      cout << "      Proceesing event " << iProgToTrain << "/" << nEvtsToTrain << "...\r" << flush;
    }

    // fill uncalibrated histograms & profiles
    hHCalFrac[0]            -> Fill(fracParVsLeadBHCal_train);
    hHCalFrac[1]            -> Fill(fracParVsSumBHCal_train);
    hECalFrac[0]            -> Fill(fracParVsLeadBEMC_train);
    hECalFrac[1]            -> Fill(fracParVsSumBEMC_train);
    hHCalDiff[0]            -> Fill(diffLeadBHCal_train);
    hHCalDiff[1]            -> Fill(diffSumBHCal_train);
    hECalDiff[0]            -> Fill(diffLeadBEMC_train);
    hECalDiff[1]            -> Fill(diffSumBEMC_train);
    hHCalEneVsPar[0]        -> Fill(ePar_train,               eLeadBHCal_train);
    pHCalEneVsPar[0]        -> Fill(ePar_train,               eLeadBHCal_train);
    hECalEneVsPar[0]        -> Fill(ePar_train,               eLeadBEMC_train);
    pECalEneVsPar[0]        -> Fill(ePar_train,               eLeadBEMC_train);
    hHCalEneVsPar[1]        -> Fill(ePar_train,               eSumBHCal_train);
    pHCalEneVsPar[1]        -> Fill(ePar_train,               eSumBHCal_train);
    hECalEneVsPar[1]        -> Fill(ePar_train,               eSumBEMC_train);
    pECalEneVsPar[1]        -> Fill(ePar_train,               eSumBEMC_train);
    hHCalFracVsPar[0]       -> Fill(ePar_train,               fracParVsLeadBHCal_train);
    pHCalFracVsPar[0]       -> Fill(ePar_train,               fracParVsLeadBHCal_train);
    hHCalFracVsPar[1]       -> Fill(ePar_train,               fracParVsSumBHCal_train);
    pHCalFracVsPar[1]       -> Fill(ePar_train,               fracParVsSumBHCal_train);
    hHCalDiffVsPar[0]       -> Fill(ePar_train,               diffLeadBHCal_train);
    pHCalDiffVsPar[0]       -> Fill(ePar_train,               diffSumBHCal_train);
    hHCalDiffVsPar[1]       -> Fill(ePar_train,               diffLeadBHCal_train);
    pHCalDiffVsPar[1]       -> Fill(ePar_train,               diffSumBHCal_train);
    hECalFracVsPar[0]       -> Fill(ePar_train,               fracParVsLeadBEMC_train);
    pECalFracVsPar[0]       -> Fill(ePar_train,               fracParVsLeadBEMC_train);
    hECalFracVsPar[1]       -> Fill(ePar_train,               fracParVsSumBEMC_train);
    pECalFracVsPar[1]       -> Fill(ePar_train,               fracParVsSumBEMC_train);
    hECalDiffVsPar[0]       -> Fill(ePar_train,               diffLeadBEMC_train);
    pECalDiffVsPar[0]       -> Fill(ePar_train,               diffSumBEMC_train);
    hECalDiffVsPar[1]       -> Fill(ePar_train,               diffLeadBEMC_train);
    pECalDiffVsPar[1]       -> Fill(ePar_train,               diffSumBEMC_train);
    hHCalVsECalFrac[0]      -> Fill(fracParVsLeadBEMC_train,  fracParVsLeadBHCal_train);
    pHCalVsECalFrac[0]      -> Fill(fracParVsLeadBEMC_train,  fracParVsLeadBHCal_train);
    hHCalVsECalFrac[1]      -> Fill(fracParVsSumBEMC_train,   fracParVsSumBHCal_train);
    pHCalVsECalFrac[1]      -> Fill(fracParVsSumBEMC_train,   fracParVsSumBHCal_train);
    hHCalVsECalDiff[0]      -> Fill(diffLeadBEMC_train,       diffLeadBHCal_train);
    pHCalVsECalDiff[0]      -> Fill(diffLeadBEMC_train,       diffLeadBHCal_train);
    hHCalVsECalDiff[1]      -> Fill(diffSumBEMC_train,        diffSumBHCal_train);
    pHCalVsECalDiff[1]      -> Fill(diffSumBEMC_train,        diffSumBHCal_train);
    hHCalFracVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC_train, fracParVsLeadBHCal_train);
    pHCalFracVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC_train, fracParVsLeadBHCal_train);
    hHCalFracVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC_train, fracParVsSumBHCal_train);
    pHCalFracVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC_train, fracParVsSumBHCal_train);
    hHCalDiffVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC_train, diffLeadBHCal_train);
    pHCalDiffVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC_train, diffLeadBHCal_train);
    hHCalDiffVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC_train, diffSumBHCal_train);
    pHCalDiffVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC_train, diffSumBHCal_train);
    hECalFracVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC_train, fracParVsLeadBEMC_train);
    pECalFracVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC_train, fracParVsLeadBEMC_train);
    hECalFracVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC_train, fracParVsSumBEMC_train);
    pECalFracVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC_train, fracParVsSumBEMC_train);
    hECalDiffVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC_train, diffLeadBEMC_train);
    pECalDiffVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC_train, diffLeadBEMC_train);
    hECalDiffVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC_train, diffSumBEMC_train);
    pECalDiffVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC_train, diffSumBEMC_train);

    // fill resolution histograms
    for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {
      const bool isInEneParBin = ((ePar_train > eneParMin[iEneBin]) && (ePar_train < eneParMax[iEneBin]));
      if (isInEneParBin) {
        hHCalEneBin[iEneBin]  -> Fill(eLeadBHCal_train);
        hHCalDiffBin[iEneBin] -> Fill(diffLeadBHCal_train);
      }
    }
  }  // end uncalibrated event loop
  cout << "    Finished uncalibrated event loop." << endl;

  // resolution calculation
  TF1      *fFitEneBin[NEneBins];
  TF1      *fFitDiffBin[NEneBins];
  double  binSigmaEne[NEneBins];
  double  valSigmaEne[NEneBins];
  double  valSigmaEneHist[NEneBins];
  double  valSigmaDiff[NEneBins];
  double  valSigmaDiffHist[NEneBins];
  double  errSigmaEne[NEneBins];
  double  errSigmaEneHist[NEneBins];
  double  errSigmaDiff[NEneBins];
  double  errSigmaDiffHist[NEneBins];
  for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {

    // normalize hisotgrams
    const double intEneBin  = hHCalEneBin[iEneBin]  -> Integral();
    const double intDiffBin = hHCalDiffBin[iEneBin] -> Integral();
    if (intEneBin > 0.)  hHCalEneBin[iEneBin] -> Scale(1. / intEneBin);
    if (intDiffBin > 0.) hHCalDiffBin[iEneBin] -> Scale(1. / intDiffBin);

    // initialize functions
    fFitEneBin[iEneBin]  = new TF1(sFitEne[iEneBin].data(),  "gaus(0)", xFitEneMin[iEneBin],  xFitEneMax[iEneBin]);
    fFitDiffBin[iEneBin] = new TF1(sFitDiff[iEneBin].data(), "gaus(0)", xFitDiffMin[iEneBin], xFitDiffMax[iEneBin]);
    fFitEneBin[iEneBin]  -> SetParameter(0, ampEneGuess[iEneBin]);
    fFitEneBin[iEneBin]  -> SetParameter(1, muEneGuess[iEneBin]);
    fFitEneBin[iEneBin]  -> SetParameter(2, sigEneGuess[iEneBin]);
    fFitDiffBin[iEneBin] -> SetParameter(0, ampDiffGuess[iEneBin]);
    fFitDiffBin[iEneBin] -> SetParameter(1, muDiffGuess[iEneBin]);
    fFitDiffBin[iEneBin] -> SetParameter(2, sigDiffGuess[iEneBin]);

    // fit histograms
    hHCalEneBin[iEneBin]  -> Fit(sFitEne[iEneBin].data(), "r");
    hHCalDiffBin[iEneBin] -> Fit(sFitDiff[iEneBin].data(), "r");

    // grab resolutions and uncertainties
    const double muEne      = fFitEneBin[iEneBin]  -> GetParameter(1);
    const double muDiff     = fFitDiffBin[iEneBin] -> GetParameter(1);
    const double sigmaEne   = fFitEneBin[iEneBin]  -> GetParameter(2);
    const double sigmaDiff  = fFitDiffBin[iEneBin] -> GetParameter(2);
    const double errMuEne   = fFitEneBin[iEneBin]  -> GetParError(1);
    const double errMuDiff  = fFitDiffBin[iEneBin] -> GetParError(1);
    const double errSigEne  = fFitEneBin[iEneBin]  -> GetParError(2);
    const double errSigDiff = fFitDiffBin[iEneBin] -> GetParError(2);
    const double perMuEne   = errMuEne / muEne;
    const double perMuDiff  = errMuDiff / muDiff;
    const double perSigEne  = errSigEne / sigmaEne;
    const double perSigDiff = errSigDiff / sigmaDiff;

    const double muHistEne      = hHCalEneBin[iEneBin]  -> GetMean();
    const double muHistDiff     = hHCalDiffBin[iEneBin] -> GetMean();
    const double sigmaHistEne   = hHCalEneBin[iEneBin]  -> GetRMS();
    const double sigmaHistDiff  = hHCalDiffBin[iEneBin] -> GetRMS();
    const double errMuHistEne   = hHCalEneBin[iEneBin]  -> GetMeanError();
    const double errMuHistDiff  = hHCalDiffBin[iEneBin] -> GetMeanError();
    const double errSigHistEne  = hHCalEneBin[iEneBin]  -> GetRMSError();
    const double errSigHistDiff = hHCalDiffBin[iEneBin] -> GetRMSError();
    const double perMuHistEne   = errMuHistEne / muHistEne;
    const double perMuHistDiff  = errMuHistDiff / muHistDiff;
    const double perSigHistEne  = errSigHistEne / sigmaHistEne;
    const double perSigHistDiff = errSigHistDiff / sigmaHistDiff;

    binSigmaEne[iEneBin]  = (eneParMin[iEneBin] - eneParMax[iEneBin]) / 2.;
    valSigmaEne[iEneBin]  = sigmaEne / muEne;
    valSigmaDiff[iEneBin] = sigmaDiff / muDiff;
    errSigmaEne[iEneBin]  = valSigmaEne[iEneBin] * TMath::Sqrt((perMuEne * perMuEne) + (perSigEne * perSigEne));
    errSigmaDiff[iEneBin] = valSigmaDiff[iEneBin] * TMath::Sqrt((perMuDiff * perMuDiff) + (perSigDiff * perSigDiff));

    valSigmaEneHist[iEneBin]  = sigmaHistEne / muHistEne;
    valSigmaDiffHist[iEneBin] = sigmaHistDiff / muHistDiff;
    errSigmaEneHist[iEneBin]  = valSigmaEneHist[iEneBin] * TMath::Sqrt((perMuHistEne * perMuHistEne) + (perSigHistEne * perSigHistEne));
    errSigmaDiffHist[iEneBin] = valSigmaDiffHist[iEneBin] * TMath::Sqrt((perMuHistDiff * perMuHistDiff) + (perSigHistDiff * perSigHistDiff));
  }
  cout << "    Normalized and fit resolution histograms." << endl;

  // create resolution graphs
  TGraphErrors *grResoEne      = new TGraphErrors(NEneBins, enePar, valSigmaEne,      binSigmaEne, errSigmaEne);
  TGraphErrors *grResoDiff     = new TGraphErrors(NEneBins, enePar, valSigmaDiff,     binSigmaEne, errSigmaDiff);
  TGraphErrors *grResoEneHist  = new TGraphErrors(NEneBins, enePar, valSigmaEneHist,  binSigmaEne, errSigmaEneHist);
  TGraphErrors *grResoDiffHist = new TGraphErrors(NEneBins, enePar, valSigmaDiffHist, binSigmaEne, errSigmaDiffHist);
  grResoEne      -> SetName("grResoEne");
  grResoDiff     -> SetName("grResoDiff");
  grResoEneHist  -> SetName("grResoEneHist");
  grResoDiffHist -> SetName("grResoDiffHist");
  cout << "    Made uncalibrated resolution graphs." << endl;

  // train tmva ---------------------------------------------------------------

  // instantiate tmva library
  Factory    *factory;
  DataLoader *loader;
  Tools::Instance();
  cout << "    Beginning calibration:" << endl;

  // create tmva factory & load data
  factory = new Factory("TMVARegression", fOutput, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");
  loader  = new DataLoader(sLoader.data());
  cout << "      Created factory and loaded data..." << endl;

  // set variables and target
  if (addSpectators) {
    for (UInt_t iSpectator = 0; iSpectator < NTmvaSpec; iSpectator++) {
      loader -> AddSpectator(sTmvaSpec[iSpectator]);
    }
  }
  for (UInt_t iVariable = 0; iVariable < NTmvaVar; iVariable++) {
    loader -> AddVariable(sTmvaVar[iVariable]);
  }
  loader -> AddTarget(sTarget);
  cout << "      Set spectators, variables, and target..." << endl;

  // add tree & prepare for training
  loader -> AddRegressionTree(ntToTrain, treeWeight);
  loader -> PrepareTrainingAndTestTree(trainCut, "nTrain_Regression=1000:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V");
  cout << "      Added tree and prepared for training..." << endl;

  // book methods
  factory -> BookMethod(loader, Types::kLD,  "LD");
  factory -> BookMethod(loader, Types::kMLP, "MLP");
  factory -> BookMethod(loader, Types::kBDT, "BDTG");
  cout << "      Booked methods..." << endl;

  // train, test, & evaluate
  factory -> TrainAllMethods();
  factory -> TestAllMethods();
  factory -> EvaluateAllMethods();
  cout << "      Trained TMVA.\n"
       << "    Finished calibration!"
       << endl;

  // apply model --------------------------------------------------------------

  // default methods to be trained + tested
  map<string, int> Use;
  for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
    const string sToUse(sMethods[iMethod].data());
    Use[sToUse] = 1;
  }
  cout << "\n==> Start TMVARegressionApplication" << endl;

  Reader *reader = new Reader( "!Color:!Silent" );
  reader -> AddVariable("eLeadBHCal",       &eLeadBHCal_apply);
  reader -> AddVariable("eLeadBEMC",        &eLeadBEMC_apply);
  reader -> AddVariable("hLeadBHCal",       &hLeadBHCal_apply);
  reader -> AddVariable("hLeadBEMC",        &hLeadBEMC_apply);
  reader -> AddVariable("fLeadBHCal",       &fLeadBHCal_apply);
  reader -> AddVariable("fLeadBEMC",        &fLeadBEMC_apply);
  reader -> AddVariable("nHitsLeadBHCal",   &nHitsLeadBHCal_apply);
  reader -> AddVariable("nHitsLeadBEMC",    &nHitsLeadBEMC_apply);
  reader -> AddVariable("eSumImage",        &eSumImage_apply);
  reader -> AddVariable("eSumSciFi",        &eSumSciFi_apply);
  reader -> AddVariable("eSumSciFiLayer1",  &eSumSciFiLayer1_apply);
  reader -> AddVariable("eSumSciFiLayer2",  &eSumSciFiLayer2_apply);
  reader -> AddVariable("eSumSciFiLayer3",  &eSumSciFiLayer3_apply);
  reader -> AddVariable("eSumSciFiLayer4",  &eSumSciFiLayer4_apply);
  reader -> AddVariable("eSumSciFiLayer5",  &eSumSciFiLayer5_apply);
  reader -> AddVariable("eSumSciFiLayer6",  &eSumSciFiLayer6_apply);
  reader -> AddVariable("eSumSciFiLayer7",  &eSumSciFiLayer7_apply);
  reader -> AddVariable("eSumSciFiLayer8",  &eSumSciFiLayer8_apply);
  reader -> AddVariable("eSumSciFiLayer9",  &eSumSciFiLayer9_apply);
  reader -> AddVariable("eSumSciFiLayer10", &eSumSciFiLayer10_apply);
  reader -> AddVariable("eSumSciFiLayer11", &eSumSciFiLayer11_apply);
  reader -> AddVariable("eSumSciFiLayer12", &eSumSciFiLayer12_apply);
  reader -> AddVariable("eSumImageLayer1",  &eSumImageLayer1_apply);
  //reader -> AddVariable("eSumImageLayer2",  &eSumImageLayer2_apply);
  reader -> AddVariable("eSumImageLayer3",  &eSumImageLayer3_apply);
  reader -> AddVariable("eSumImageLayer4",  &eSumImageLayer4_apply);
  //reader -> AddVariable("eSumImageLayer5",  &eSumImageLayer5_apply);
  reader -> AddVariable("eSumImageLayer6",  &eSumImageLayer6_apply);

  // book method(s)
  for (map<string, int>::iterator itMethod = Use.begin(); itMethod != Use.end(); itMethod++) {
    if (itMethod -> second) {
      string methodName = string(itMethod -> first) + " method";
      string weightfile = sLoader + "/weights/" + STmvaPrefix + "_" + string(itMethod -> first) + ".weights.xml";
      reader->BookMVA(methodName, weightfile);
    }
  }  // end method loop

  // for tmva histogram binning
  const UInt_t  nTmvaBins(100);
  const float rTmvaBins[NRange] = {-100., 600.};

  // Book tmva histograms
  Int_t  nTmvaHist(-1);
  TH1   *hTMVA[NTmvaHistMax];
  for (map<string, int>::iterator itMethod = Use.begin(); itMethod != Use.end(); itMethod++) {
    string  sName  = string(itMethod -> first.c_str());
    string  sTitle = string(itMethod -> first) + " method";
    TH1     *hNew   = new TH1F(sName.data(), sTitle.data(), nTmvaBins, rTmvaBins[0], rTmvaBins[1]);
    if (!hNew) {
      cerr << "PANIC: couldn't create TMVA histogram #" << nTmvaHist << "! Aborting code execution!\n" << endl;
      return;
    } else {
      if (itMethod -> second) hTMVA[++nTmvaHist] = hNew;
    }
  }  // end method loop
  nTmvaHist++;

  // get number of events for application
  const Long64_t nEvtsToApply = ntToApply -> GetEntries();

  // begin application event loop
  TStopwatch stopwatch;
  cout << "--- Processing: " << nEvtsToApply << " events" << endl;

  nBytes = 0;
  stopwatch.Start();
  for (Long64_t iEvtToApply = 0; iEvtToApply < nEvtsToApply; iEvtToApply++) {

    // announce progress
    if (iEvtToApply % 1000 == 0) {
      cout << "--- ... Processing event: " << iEvtToApply << endl;
    }

    const Long64_t bytes = ntToApply -> GetEntry(iEvtToApply);
    if (bytes < 0.) {
      cerr << "WARNING something wrong with event " << iEvtToApply << "! Aborting loop!" << endl;
      break;
    }
    nBytes += bytes;

    // loop over methods
    for (Int_t iTmvaHist = 0; iTmvaHist < nTmvaHist; iTmvaHist++) {

      // grab regression target
      TString title  = hTMVA[iTmvaHist] -> GetTitle();
      float target = (reader -> EvaluateRegression(title))[0];
      hTMVA[iTmvaHist] -> Fill(target);

      // check for method
      Int_t method = -1;
      for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
        bool isMethod = title.Contains(sMethods[iMethod].data());
        if (isMethod) {
          method = iMethod;
          break;
        }
      }  // end method loop

      // check for ecal energy
      const bool methodExists     = (method > -1);
      const bool isInECalEneRange = ((eLeadBEMC_apply > eneECalRange[0]) && (eLeadBEMC_apply < eneECalRange[1]));
      if (doECalCut && !isInECalEneRange) continue;

      // fill resolution histograms
      if (methodExists) {
        for (UInt_t iCalibBin = 0; iCalibBin < NCalibBins; iCalibBin++) {
          const bool isInEneCalibBin = ((ePar_apply > eneCalibMin[iCalibBin]) && (ePar_apply < eneCalibMax[iCalibBin]));
          if (isInEneCalibBin) {
            hHCalCalibBin[method][iCalibBin] -> Fill(target);
          }
        }  // end energy bin loop
        hCalibCalibVsPar[method]  -> Fill(ePar_apply,      target);
        hHCalCalibVsPar[method]   -> Fill(ePar_apply,      eLeadBHCal_apply);
        hHCalCalibVsCalib[method] -> Fill(target,          eLeadBHCal_apply);
        hHCalCalibVsECal[method]  -> Fill(eLeadBEMC_apply, eLeadBHCal_apply);
        hECalCalibVsPar[method]   -> Fill(ePar_apply,      eLeadBEMC_apply);
        hECalCalibVsCalib[method] -> Fill(target,          eLeadBEMC_apply);
      }  // end if (methodExists)
    }  // end method loop
  }  // end event loop
  stopwatch.Stop();

  // announce end of event loop
  cout << "--- End of event loop: ";
  cout << "\n    Application finished!" << endl;
  stopwatch.Print();

  // calculate calibrated resolution ------------------------------------------

  // for graphs
  double binSigmaCalib[NMethods][NCalibBins];
  double valMuCalib[NMethods][NCalibBins];
  double valMuCalibHist[NMethods][NCalibBins];
  double valSigmaCalib[NMethods][NCalibBins];
  double valSigmaCalibHist[NMethods][NCalibBins];
  double errMuCalib[NMethods][NCalibBins];
  double errMuCalibHist[NMethods][NCalibBins];
  double errSigmaCalib[NMethods][NCalibBins];
  double errSigmaCalibHist[NMethods][NCalibBins];

  // resolution calculation
  TF1          *fFitCalibBin[NMethods][NCalibBins];
  TGraphErrors *grLineCalib[NMethods];
  TGraphErrors *grLineCalibHist[NMethods];
  TGraphErrors *grResoCalib[NMethods];
  TGraphErrors *grResoCalibHist[NMethods];
  for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
    for (UInt_t iCalibBin = 0; iCalibBin < NCalibBins; iCalibBin++) {

      // normalize hisotgrams
      const double intCalibBin = hHCalCalibBin[iMethod][iCalibBin]  -> Integral();
      if (intCalibBin > 0.) hHCalCalibBin[iMethod][iCalibBin] -> Scale(1. / intCalibBin);

      // make name
      TString sFitCalib(sFitCalibBase[iCalibBin].data());
      sFitCalib.Append("_");
      sFitCalib.Append(sMethods[iMethod]);

      // initialize functions
      fFitCalibBin[iMethod][iCalibBin] = new TF1(sFitCalib.Data(),  "gaus(0)", xFitCalibMin[iCalibBin],  xFitCalibMax[iCalibBin]);
      fFitCalibBin[iMethod][iCalibBin] -> SetParameter(0, ampCalibGuess[iCalibBin]);
      fFitCalibBin[iMethod][iCalibBin] -> SetParameter(1, muCalibGuess[iCalibBin]);
      fFitCalibBin[iMethod][iCalibBin] -> SetParameter(2, sigCalibGuess[iCalibBin]);

      // fit histograms
      hHCalCalibBin[iMethod][iCalibBin] -> Fit(sFitCalib.Data(), "rQ");

      // grab resolutions and uncertainties
      const double mu     = fFitCalibBin[iMethod][iCalibBin]  -> GetParameter(1);
      const double sigma  = fFitCalibBin[iMethod][iCalibBin]  -> GetParameter(2);
      const double errMu  = fFitCalibBin[iMethod][iCalibBin]  -> GetParError(1);
      const double errSig = fFitCalibBin[iMethod][iCalibBin]  -> GetParError(2);
      const double perMu  = errMu / mu;
      const double perSig = errSig / sigma;

      const double muHist     = hHCalCalibBin[iMethod][iCalibBin]  -> GetMean();
      const double sigmaHist  = hHCalCalibBin[iMethod][iCalibBin]  -> GetRMS();
      const double errMuHist  = hHCalCalibBin[iMethod][iCalibBin]  -> GetMeanError();
      const double errSigHist = hHCalCalibBin[iMethod][iCalibBin]  -> GetRMSError();
      const double perMuHist  = errMuHist / muHist;
      const double perSigHist = errSigHist / sigmaHist;

      // set fit values
      binSigmaCalib[iMethod][iCalibBin] = (eneParMin[iCalibBin] - eneParMax[iCalibBin]) / 2.;
      valMuCalib[iMethod][iCalibBin]    = mu;
      valSigmaCalib[iMethod][iCalibBin] = sigma / mu;
      errMuCalib[iMethod][iCalibBin]    = errMu;
      errSigmaCalib[iMethod][iCalibBin] = valSigmaCalib[iMethod][iCalibBin] * TMath::Sqrt((perMu * perMu) + (perSig * perSig));

      // set histogram values
      valMuCalibHist[iMethod][iCalibBin]    = muHist;
      valSigmaCalibHist[iMethod][iCalibBin] = sigmaHist / muHist;
      errMuCalibHist[iMethod][iCalibBin]    = errMuHist;
      errSigmaCalibHist[iMethod][iCalibBin] = valSigmaCalibHist[iMethod][iCalibBin] * TMath::Sqrt((perMuHist * perMuHist) + (perSigHist * perSigHist));
    }
    cout << "    Fit resolution histograms." << endl;

    // make name
    TString sGraphLineCalib("grLineCalib");
    TString sGraphLineCalibHist("grLineCalibHist");
    TString sGraphResoCalib("grResoCalib");
    TString sGraphResoCalibHist("grResoCalibHist");
    sGraphLineCalib.Append("_");
    sGraphLineCalibHist.Append("_");
    sGraphResoCalib.Append("_");
    sGraphResoCalibHist.Append("_");
    sGraphLineCalib.Append(sMethods[iMethod].data());
    sGraphLineCalibHist.Append(sMethods[iMethod].data());
    sGraphResoCalib.Append(sMethods[iMethod].data());
    sGraphResoCalibHist.Append(sMethods[iMethod].data());

    // create resolution graphs
    grLineCalib[iMethod]     = new TGraphErrors(NCalibBins, enePar, valMuCalib[iMethod],        binSigmaCalib[iMethod], errMuCalib[iMethod]);
    grLineCalibHist[iMethod] = new TGraphErrors(NCalibBins, enePar, valMuCalibHist[iMethod],    binSigmaCalib[iMethod], errMuCalibHist[iMethod]);
    grResoCalib[iMethod]     = new TGraphErrors(NCalibBins, enePar, valSigmaCalib[iMethod],     binSigmaCalib[iMethod], errSigmaCalib[iMethod]);
    grResoCalibHist[iMethod] = new TGraphErrors(NCalibBins, enePar, valSigmaCalibHist[iMethod], binSigmaCalib[iMethod], errSigmaCalibHist[iMethod]);
    grLineCalib[iMethod]     -> SetName(sGraphLineCalib.Data());
    grLineCalibHist[iMethod] -> SetName(sGraphLineCalibHist.Data());
    grResoCalib[iMethod]     -> SetName(sGraphResoCalib.Data());
    grResoCalibHist[iMethod] -> SetName(sGraphResoCalibHist.Data());
  }  // end method loop

  // save output and close ----------------------------------------------------

  // save histograms
  TDirectory *dUncal = (TDirectory*) fOutput -> mkdir("uncalibrated");
  TDirectory *dCalib = (TDirectory*) fOutput -> mkdir("calibrated");
  TDirectory *dReso  = (TDirectory*) fOutput -> mkdir("resolution");
  TDirectory *dTmva = (TDirectory*) fOutput -> mkdir("tmva");
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    if (!isCalibrated[iHist]) {
      dUncal -> cd();
    }
    hHCalFrac[iHist]            -> Write();
    hHCalDiff[iHist]            -> Write();
    hECalFrac[iHist]            -> Write();
    hECalDiff[iHist]            -> Write();
    hHCalEneVsPar[iHist]        -> Write();
    pHCalEneVsPar[iHist]        -> Write();
    hECalEneVsPar[iHist]        -> Write();
    pECalEneVsPar[iHist]        -> Write();
    hHCalFracVsPar[iHist]       -> Write();
    pHCalFracVsPar[iHist]       -> Write();
    hHCalDiffVsPar[iHist]       -> Write();
    pHCalDiffVsPar[iHist]       -> Write();
    hECalFracVsPar[iHist]       -> Write();
    pECalFracVsPar[iHist]       -> Write();
    hECalDiffVsPar[iHist]       -> Write();
    pECalDiffVsPar[iHist]       -> Write();
    hHCalVsECalFrac[iHist]      -> Write();
    pHCalVsECalFrac[iHist]      -> Write();
    hHCalVsECalDiff[iHist]      -> Write();
    pHCalVsECalDiff[iHist]      -> Write();
    hHCalFracVsTotalFrac[iHist] -> Write();
    pHCalFracVsTotalFrac[iHist] -> Write();
    hHCalDiffVsTotalFrac[iHist] -> Write();
    pHCalDiffVsTotalFrac[iHist] -> Write();
    hECalFracVsTotalFrac[iHist] -> Write();
    pECalFracVsTotalFrac[iHist] -> Write();
    hECalDiffVsTotalFrac[iHist] -> Write();
    pECalDiffVsTotalFrac[iHist] -> Write();
  }

  dReso          -> cd();
  grResoEne      -> Write();
  grResoDiff     -> Write();
  grResoEneHist  -> Write();
  grResoDiffHist -> Write();
  for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {
    hHCalEneBin[iEneBin]  -> Write();
    hHCalDiffBin[iEneBin] -> Write();
    fFitEneBin[iEneBin]   -> Write();
    fFitDiffBin[iEneBin]  -> Write();
  }

  dReso -> cd();
  for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
    hCalibCalibVsPar[iMethod]  -> Write();
    hHCalCalibVsPar[iMethod]   -> Write();
    hHCalCalibVsCalib[iMethod] -> Write();
    hHCalCalibVsECal[iMethod]  -> Write();
    hECalCalibVsPar[iMethod]   -> Write();
    hECalCalibVsCalib[iMethod] -> Write();
    grLineCalib[iMethod]       -> Write();
    grLineCalibHist[iMethod]   -> Write();
    grResoCalib[iMethod]       -> Write();
    grResoCalibHist[iMethod]   -> Write();
    for (UInt_t iCalibBin = 0; iCalibBin < NCalibBins; iCalibBin++) {
      hHCalCalibBin[iMethod][iCalibBin] -> Write();
      fFitCalibBin[iMethod][iCalibBin]  -> Write();
    }
  }  // end method loop
  cout << "    Saved histograms." << endl;

  // close files
  fOutput  -> cd();
  fOutput  -> Close();
  fInTrain -> cd();
  fInTrain -> Close();
  fInApply -> cd();
  fInApply -> Close();
  cout << "  Finished BHCal calibration script!\n" << endl;

  // delete tmva objects and exit
  delete factory;
  delete loader;
  delete reader;
  return;

}

// end ------------------------------------------------------------------------
