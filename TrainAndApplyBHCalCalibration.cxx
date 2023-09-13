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
// dataframe related classes
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>

// make common namespaces implicit
using namespace std;
using namespace TMVA;
using namespace ROOT;
using namespace ROOT::VecOps;

// global constants
static const UInt_t NTxt(3);
static const UInt_t NVtx(4);
static const UInt_t NHist(4);
static const UInt_t NRange(2);
static const UInt_t NEneBins(10);
static const UInt_t NTmvaVar(28);
static const UInt_t NTmvaSpec(1);

// default arguments
static const string SInDef("../performance/eicrecon_output/single_particles/merged/forPerformanceStudy.withIndividualECalLayers_includedEPar7.e110th45n20Kneu.d20m7y2023.plugin.root");
static const string SOutDef("StreamlineTest_Change0_PrunedTraining.train.root");
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
  const string sLoader("StreamlineTest_Baseline");
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
    "eSumImageLayer2",
    "eSumImageLayer3",
    "eSumImageLayer4",
    "eSumImageLayer5",
    "eSumImageLayer6"
  };
  const string sTmvaSpec[NTmvaSpec] = {""};

  // histogram parameters
  const bool    isCalibrated[NHist] = {false, false, true, true};
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

  // generic resolution parameters
  const double enePar[NEneBins]    = {2.,  3.,  4.,  5.,  6.,  8,   10.,  12.,  16.,  20.};
  const double eneParMin[NEneBins] = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 9.5,  11.5, 13.5, 18.5};
  const double eneParMax[NEneBins] = {2.5, 3.5, 4.5, 5.5, 6.5, 9.5, 11.5, 13.5, 18.5, 21.5};

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

  // load input ---------------------------------------------------------------

  // open files
  TFile *fInput  = new TFile(sInput.data(), "read");
  TFile *fOutput = new TFile(sOutput.data(), "recreate");
  if (!fInput || !fOutput) {
    cerr << "PANIC: couldn't open a file!\n"
         << "        fInput = " << fInput << ", fOutput = " << fOutput << "\n"
         << endl;
    return;
  }
  cout << "    Opened files:\n"
       << "      fInput  = " << sInput.data() << "\n"
       << "      fOutput = " << sOutput.data()
       << endl;

  // grab input tuple
  TNtuple *ntToCalibrate = (TNtuple*) fInput -> Get(sTuple.data());
  if (!ntToCalibrate) {
    cerr << "PANIC: couldn't grab input tuple!\n"
         << "       name  = " << sTuple        << "\n"
         << "       tuple = " << ntToCalibrate << "\n"
         << endl;
    return;
  }
  cout << "    Grabbed input tuple:\n"
       << "      tuple = " << sTuple
       << endl;

  // declare tuple leaves
  float ePar;
  float fracParVsLeadBHCal;
  float fracParVsLeadBEMC;
  float fracParVsSumBHCal;
  float fracParVsSumBEMC;
  float fracLeadBHCalVsBEMC;
  float fracSumBHCalVsBEMC;
  float eLeadBHCal;
  float eLeadBEMC;
  float eSumBHCal;
  float eSumBEMC;
  float diffLeadBHCal;
  float diffLeadBEMC;
  float diffSumBHCal;
  float diffSumBEMC;
  float nHitsLeadBHCal;
  float nHitsLeadBEMC;
  float nClustBHCal;
  float nClustBEMC;
  float hLeadBHCal;
  float hLeadBEMC;
  float fLeadBHCal;
  float fLeadBEMC;
  float eLeadImage;
  float eSumImage;
  float eLeadSciFi;
  float eSumSciFi;
  float nClustImage;
  float nClustSciFi;
  float hLeadImage;
  float hLeadSciFi;
  float fLeadImage;
  float fLeadSciFi;
  float eSumSciFiLayer1;
  float eSumSciFiLayer2;
  float eSumSciFiLayer3;
  float eSumSciFiLayer4;
  float eSumSciFiLayer5;
  float eSumSciFiLayer6;
  float eSumSciFiLayer7;
  float eSumSciFiLayer8;
  float eSumSciFiLayer9;
  float eSumSciFiLayer10;
  float eSumSciFiLayer11;
  float eSumSciFiLayer12;
  float eSumImageLayer1;
  float eSumImageLayer2;
  float eSumImageLayer3;
  float eSumImageLayer4;
  float eSumImageLayer5;
  float eSumImageLayer6;

  // set tuple branches
  ntToCalibrate -> SetBranchAddress("ePar",                &ePar);
  ntToCalibrate -> SetBranchAddress("fracParVsLeadBHCal",  &fracParVsLeadBHCal);
  ntToCalibrate -> SetBranchAddress("fracParVsLeadBEMC",   &fracParVsLeadBEMC);
  ntToCalibrate -> SetBranchAddress("fracParVsSumBHCal",   &fracParVsSumBHCal);
  ntToCalibrate -> SetBranchAddress("fracParVsSumBEMC",    &fracParVsSumBEMC);
  ntToCalibrate -> SetBranchAddress("fracLeadBHCalVsBEMC", &fracLeadBHCalVsBEMC);
  ntToCalibrate -> SetBranchAddress("fracSumBHCalVsBEMC",  &fracSumBHCalVsBEMC);
  ntToCalibrate -> SetBranchAddress("eLeadBHCal",          &eLeadBHCal);
  ntToCalibrate -> SetBranchAddress("eLeadBEMC",           &eLeadBEMC);
  ntToCalibrate -> SetBranchAddress("eSumBHCal",           &eSumBHCal);
  ntToCalibrate -> SetBranchAddress("eSumBEMC",            &eSumBEMC);
  ntToCalibrate -> SetBranchAddress("diffLeadBHCal",       &diffLeadBHCal);
  ntToCalibrate -> SetBranchAddress("diffLeadBEMC",        &diffLeadBEMC);
  ntToCalibrate -> SetBranchAddress("diffSumBHCal",        &diffSumBHCal);
  ntToCalibrate -> SetBranchAddress("diffSumBEMC",         &diffSumBEMC);
  ntToCalibrate -> SetBranchAddress("nHitsLeadBHCal",      &nHitsLeadBHCal);
  ntToCalibrate -> SetBranchAddress("nHitsLeadBEMC",       &nHitsLeadBEMC);
  ntToCalibrate -> SetBranchAddress("nClustBHCal",         &nClustBHCal);
  ntToCalibrate -> SetBranchAddress("nClustBEMC",          &nClustBEMC);
  ntToCalibrate -> SetBranchAddress("hLeadBHCal",          &hLeadBHCal);
  ntToCalibrate -> SetBranchAddress("hLeadBEMC",           &hLeadBEMC);
  ntToCalibrate -> SetBranchAddress("fLeadBHCal",          &fLeadBHCal);
  ntToCalibrate -> SetBranchAddress("fLeadBEMC",           &fLeadBEMC);
  ntToCalibrate -> SetBranchAddress("eLeadImage",          &eLeadImage);
  ntToCalibrate -> SetBranchAddress("eSumImage",           &eSumImage);
  ntToCalibrate -> SetBranchAddress("eLeadSciFi",          &eLeadSciFi);
  ntToCalibrate -> SetBranchAddress("eSumSciFi",           &eSumSciFi);
  ntToCalibrate -> SetBranchAddress("nClustImage",         &nClustImage);
  ntToCalibrate -> SetBranchAddress("nClustSciFi",         &nClustSciFi);
  ntToCalibrate -> SetBranchAddress("hLeadImage",          &hLeadImage);
  ntToCalibrate -> SetBranchAddress("hLeadSciFi",          &hLeadSciFi);
  ntToCalibrate -> SetBranchAddress("fLeadImage",          &fLeadImage);
  ntToCalibrate -> SetBranchAddress("fLeadSciFi",          &fLeadSciFi);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer1",     &eSumSciFiLayer1);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer2",     &eSumSciFiLayer2);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer3",     &eSumSciFiLayer3);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer4",     &eSumSciFiLayer4);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer5",     &eSumSciFiLayer5);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer6",     &eSumSciFiLayer6);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer7",     &eSumSciFiLayer7);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer8",     &eSumSciFiLayer8);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer9",     &eSumSciFiLayer9);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer10",    &eSumSciFiLayer10);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer11",    &eSumSciFiLayer11);
  ntToCalibrate -> SetBranchAddress("eSumSciFiLayer12",    &eSumSciFiLayer12);
  ntToCalibrate -> SetBranchAddress("eSumImageLayer1",     &eSumImageLayer1);
  ntToCalibrate -> SetBranchAddress("eSumImageLayer2",     &eSumImageLayer2);
  ntToCalibrate -> SetBranchAddress("eSumImageLayer3",     &eSumImageLayer3);
  ntToCalibrate -> SetBranchAddress("eSumImageLayer4",     &eSumImageLayer4);
  ntToCalibrate -> SetBranchAddress("eSumImageLayer5",     &eSumImageLayer5);
  ntToCalibrate -> SetBranchAddress("eSumImageLayer6",     &eSumImageLayer6);
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
  const Ssiz_t  nEneBins(41);
  const Ssiz_t  nDiffBins(700);
  const Ssiz_t  nFracBins(305);
  const float rEneBins[NRange]  = {-1.,   40.};
  const float rDiffBins[NRange] = {-1.5,  5.5};
  const float rFracBins[NRange] = {-0.05, 3.};

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
  cout << "    Declared output histograms." << endl;

  // loop over ntuple entries -------------------------------------------------

  // prepare for uncalibrated tuple loop
  Long64_t nEvts = ntToCalibrate -> GetEntries();
  cout << "    Looping over uncalibrated tuple: " << nEvts << " events to process." << endl;

  Long64_t nBytes(0);
  for (Long64_t iEvt = 0; iEvt < nEvts; iEvt++) {

    const Long64_t bytes = ntToCalibrate -> GetEntry(iEvt);
    if (bytes < 0.) {
      cerr << "WARNING something wrong with event " << iEvt << "! Aborting loop!" << endl;
      break;
    }
    nBytes += bytes;

    // announce progress
    const Long64_t iProg = iEvt + 1;
    if (iProg == nEvts) {
      cout << "      Proceesing event " << iProg << "/" << nEvts << "..." << endl;
    } else {
      cout << "      Proceesing event " << iProg << "/" << nEvts << "...\r" << flush;
    }

    // fill uncalibrated histograms & profiles
    hHCalFrac[0]            -> Fill(fracParVsLeadBHCal);
    hHCalFrac[1]            -> Fill(fracParVsSumBHCal);
    hECalFrac[0]            -> Fill(fracParVsLeadBEMC);
    hECalFrac[1]            -> Fill(fracParVsSumBEMC);
    hHCalDiff[0]            -> Fill(diffLeadBHCal);
    hHCalDiff[1]            -> Fill(diffSumBHCal);
    hECalDiff[0]            -> Fill(diffLeadBEMC);
    hECalDiff[1]            -> Fill(diffSumBEMC);
    hHCalEneVsPar[0]        -> Fill(ePar,               eLeadBHCal);
    pHCalEneVsPar[0]        -> Fill(ePar,               eLeadBHCal);
    hECalEneVsPar[0]        -> Fill(ePar,               eLeadBEMC);
    pECalEneVsPar[0]        -> Fill(ePar,               eLeadBEMC);
    hHCalEneVsPar[1]        -> Fill(ePar,               eSumBHCal);
    pHCalEneVsPar[1]        -> Fill(ePar,               eSumBHCal);
    hECalEneVsPar[1]        -> Fill(ePar,               eSumBEMC);
    pECalEneVsPar[1]        -> Fill(ePar,               eSumBEMC);
    hHCalFracVsPar[0]       -> Fill(ePar,               fracParVsLeadBHCal);
    pHCalFracVsPar[0]       -> Fill(ePar,               fracParVsLeadBHCal);
    hHCalFracVsPar[1]       -> Fill(ePar,               fracParVsSumBHCal);
    pHCalFracVsPar[1]       -> Fill(ePar,               fracParVsSumBHCal);
    hHCalDiffVsPar[0]       -> Fill(ePar,               diffLeadBHCal);
    pHCalDiffVsPar[0]       -> Fill(ePar,               diffSumBHCal);
    hHCalDiffVsPar[1]       -> Fill(ePar,               diffLeadBHCal);
    pHCalDiffVsPar[1]       -> Fill(ePar,               diffSumBHCal);
    hECalFracVsPar[0]       -> Fill(ePar,               fracParVsLeadBEMC);
    pECalFracVsPar[0]       -> Fill(ePar,               fracParVsLeadBEMC);
    hECalFracVsPar[1]       -> Fill(ePar,               fracParVsSumBEMC);
    pECalFracVsPar[1]       -> Fill(ePar,               fracParVsSumBEMC);
    hECalDiffVsPar[0]       -> Fill(ePar,               diffLeadBEMC);
    pECalDiffVsPar[0]       -> Fill(ePar,               diffSumBEMC);
    hECalDiffVsPar[1]       -> Fill(ePar,               diffLeadBEMC);
    pECalDiffVsPar[1]       -> Fill(ePar,               diffSumBEMC);
    hHCalVsECalFrac[0]      -> Fill(fracParVsLeadBEMC,  fracParVsLeadBHCal);
    pHCalVsECalFrac[0]      -> Fill(fracParVsLeadBEMC,  fracParVsLeadBHCal);
    hHCalVsECalFrac[1]      -> Fill(fracParVsSumBEMC,   fracParVsSumBHCal);
    pHCalVsECalFrac[1]      -> Fill(fracParVsSumBEMC,   fracParVsSumBHCal);
    hHCalVsECalDiff[0]      -> Fill(diffLeadBEMC,       diffLeadBHCal);
    pHCalVsECalDiff[0]      -> Fill(diffLeadBEMC,       diffLeadBHCal);
    hHCalVsECalDiff[1]      -> Fill(diffSumBEMC,        diffSumBHCal);
    pHCalVsECalDiff[1]      -> Fill(diffSumBEMC,        diffSumBHCal);
    hHCalFracVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC, fracParVsLeadBHCal);
    pHCalFracVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC, fracParVsLeadBHCal);
    hHCalFracVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC, fracParVsSumBHCal);
    pHCalFracVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC, fracParVsSumBHCal);
    hHCalDiffVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC, diffLeadBHCal);
    pHCalDiffVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC, diffLeadBHCal);
    hHCalDiffVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC, diffSumBHCal);
    pHCalDiffVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC, diffSumBHCal);
    hECalFracVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC, fracParVsLeadBEMC);
    pECalFracVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC, fracParVsLeadBEMC);
    hECalFracVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC, fracParVsSumBEMC);
    pECalFracVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC, fracParVsSumBEMC);
    hECalDiffVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC, diffLeadBEMC);
    pECalDiffVsTotalFrac[0] -> Fill(fracSumBHCalVsBEMC, diffLeadBEMC);
    hECalDiffVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC, diffSumBEMC);
    pECalDiffVsTotalFrac[1] -> Fill(fracSumBHCalVsBEMC, diffSumBEMC);

    // fill resolution histograms
    for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {
      const bool isInEneParBin = ((ePar > eneParMin[iEneBin]) && (ePar < eneParMax[iEneBin]));
      if (isInEneParBin) {
        hHCalEneBin[iEneBin]  -> Fill(eLeadBHCal);
        hHCalDiffBin[iEneBin] -> Fill(diffLeadBHCal);
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
  cout << "    Made resolution graphs." << endl;

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
  loader -> AddRegressionTree(ntToCalibrate, treeWeight);
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

  /* TODO apply trained model here */

  // save output and close ----------------------------------------------------

  // save histograms
  TDirectory *dUncal = (TDirectory*) fOutput -> mkdir("Uncalibrated");
  TDirectory *dCalib = (TDirectory*) fOutput -> mkdir("Calibrated");
  TDirectory *dReso  = (TDirectory*) fOutput -> mkdir("Resolution");
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    if (isCalibrated[iHist]) {
      dCalib -> cd();
    } else {
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
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Finished BHCal calibration script!\n" << endl;

  // delete tmva objects and exit
  delete factory;
  delete loader;
  return;

}

// end ------------------------------------------------------------------------
