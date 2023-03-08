// ----------------------------------------------------------------------------
// 'DoHCalCalibration.C'
// Derek Anderson
// 03.02.2023
//
// Use this to train TMVA on the output of the
// JCalibrateHCal* (or PCalibrateHCal*) plugins
// and calibrate the BHCal response.
//
// Note: 'fConfig' sets which BEMC configuration
// is being used.
//   fConfig = 0: SciGlass BEMC
//   fConfig = 1: Imaging BEMC
//   fConfig = 2: Default (see "parse configuration" block)
// ----------------------------------------------------------------------------

// standard c includes
#include <map>
#include <string>
#include <cstdlib>
#include <iostream>
// root includes
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TCut.h>
#include <TFile.h>
#include <TTree.h>
#include <TError.h>
#include <TString.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TPaveText.h>
#include <TDirectory.h>
#include <TObjString.h>
#include <TGraphErrors.h>
// tmva includes
#include <TMVA/Tools.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/TMVARegGui.h>

using namespace std;

// global constants
static const UInt_t  NHist(4);
static const UInt_t  NRange(2);
static const UInt_t  NEneBins(11);
static const UInt_t  NVarSci(3);
static const UInt_t  NVarIma(3);
static const UInt_t  NSpecSci(1);
static const UInt_t  NSpecIma(1);

// default arguments
static const UInt_t  FConfigDef(2);
static const Bool_t  DoTmvaDef(false);
static const TString SInDef("eicrecon_output/merged/forECalStudy.bhcalTestBeamConfig.e1t20th35145n5KeaPim.d7m3y2023.plugin.root");
static const TString SOutDef("test.root");
static const TString STupleDef("ntForCalibration");



void DoHCalCalibration(const UInt_t fConfig = FConfigDef, const Bool_t doTMVA = DoTmvaDef, const TString sInput = SInDef, const TString sOutput = SOutDef, const TString sTuple = STupleDef) {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning BHCal calibration script..." << endl;

  // tmva parameters
  const Float_t treeWeight(1.0);
  const TString sTarget("ePar");
  const Bool_t  addSpectators(false);
  const TString sVarSci[NVarSci]   = {"eLeadBHCal", "eLeadBEMC", "eLeadBHCal+eLeadBEMC"};
  const TString sVarIma[NVarIma]   = {"eLeadBHCal", "eLeadBEMC", "eLeadBHCal+eLeadBEMC"};
  const TString sSpecSci[NSpecSci] = {"eLeadBHCal/ePar"};
  const TString sSpecIma[NSpecIma] = {"eLeadBHCal/ePar"};
  const TCut    trainCutSci("");
  const TCut    trainCutIma("");

  // histogram parameters
  const Bool_t  isCalibrated[NHist]  = {false, false, true, true};
  const UInt_t  fColEneBin[NEneBins] = {799, 809, 634, 899, 909, 618, 879, 889, 602, 859, 869};
  const UInt_t  fMarEneBin[NEneBins] = {24,  26,  32,  25,  27,  28,  30,  24,  26,  32,  25};
  const TString sHCalEne[NEneBins]   = {"hHCalEne_ene1",  "hHCalEne_ene2",   "hHCalEne_ene3",   "hHCalEne_ene4",   "hHCalEne_ene5",  "hHCalEne_ene6",
                                        "hHCalEne_ene8",  "hHCalEne_ene10",  "hHCalEne_ene12",  "hHCalEne_ene16",  "hHCalEne_ene20"};
  const TString sHCalDiff[NEneBins]  = {"hHCalDiff_ene1", "hHCalDiff_ene2",  "hHCalDiff_ene3",  "hHCalDiff_ene4",  "hHCalDiff_ene5", "hHCalDiff_ene6",
                                        "hHCalDiff_ene8", "hHCalDiff_ene10", "hHCalDiff_ene12", "hHCalDiff_ene16", "hHCalDiff_ene20"};
  const TString sEneTitleX("E_{lead}^{BHCal} [GeV]");
  const TString sDiffTitleX("#DeltaE / E_{par}");

  // generic resolution parameters
  const Double_t enePar[NEneBins]       = {1.,  2.,  3.,  4.,  5.,  6.,  8.,  10.,  12.,  16.,  20.};
  const Double_t eneParMin[NEneBins]    = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 9.5,  10.5, 13.5, 18.5};
  const Double_t eneParMax[NEneBins]    = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 9.5, 10.5, 13.5, 18.5, 21.5};

  // reco vs. par ene resolution parameters
  const Double_t xFitEneMin[NEneBins]   = {0.,  0.,  0.,  1.,  2.,  3.,  3.,  4.,   6.,   7.,   11.};
  const Double_t xFitEneMax[NEneBins]   = {2.,  2.,  4.,  5.,  6.,  7.,  9.,  12.,  14.,  17.,  21.};
  const Double_t ampEneGuess[NEneBins]  = {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,   1.,   1.,   1.};
  const Double_t muEneGuess[NEneBins]   = {0.5, 1.5, 2.,  3.,  4.,  5.,  6.,  8.,   10.,  12.,  15.};
  const Double_t sigEneGuess[NEneBins]  = {1.,  1.,  2.,  2.,  2.,  2.,  3.,  4.,   4.,   5.,   6.};
  const TString  sFitEne[NEneBins]       = {"fFitEne_ene1",  "fFitEne_ene2",   "fFitEne_ene3",   "fFitEne_ene4",   "fFitEne_ene5",  "fFitEne_ene6",
                                            "fFitEne_ene8",  "fFitEne_ene10",  "fFitEne_ene12",  "fFitEne_ene16",  "fFitEne_ene20"};

  // diff vs. par ene resolution parameters
  const Double_t xFitDiffMin[NEneBins]  = {-1, -1.,  -1., -1., -1., -1., -1., -1.,  -1.,  -1.,  -1.};
  const Double_t xFitDiffMax[NEneBins]  = {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,   1.,   1.,   1.};
  const Double_t ampDiffGuess[NEneBins] = {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,   1.,   1.,   1.};
  const Double_t muDiffGuess[NEneBins]  = {1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,   1.,   1.,   1.};
  const Double_t sigDiffGuess[NEneBins] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,  0.1,  0.1,  0.1};
  const TString  sFitDiff[NEneBins]     = {"fFitDiff_ene1", "fFitDiff_ene2",  "fFitDiff_ene3",  "fFitDiff_ene4",  "fFitDiff_ene5", "fFitDiff_ene6",
                                           "fFitDiff_ene8", "fFitDiff_ene10", "fFitDiff_ene12", "fFitDiff_ene16", "fFitDiff_ene20"};

  // style parameters
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCenter(1);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.2);

  // parse configuration
  Bool_t  inSciGlassConfig;
  TString sTupleDir;
  switch (fConfig) {
    case 0:
      inSciGlassConfig = true;
      sTupleDir        = "JCalibrateHCalWithSciGlass/";
      break;
    case 1:
      inSciGlassConfig = false;
      sTupleDir        = "JCalibrateHCalWithImaging/";
      break;
    default:
      inSciGlassConfig = true;
      sTupleDir        = "JCalibrateHCal/";
      break;
  }

  if (inSciGlassConfig) {
    cout << "    Using SciGlass configuration..." << endl;
  } else {
    cout << "    Using imaging configuration..." << endl;
  }

  // open files
  TFile *fInput  = new TFile(sInput.Data(), "read");
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fInput || !fOutput) {
    cerr << "PANIC: couldn't open a file!\n"
         << "        fInput = " << fInput << ", fOutput = " << fOutput << "\n"
         << endl;
    return;
  }
  cout << "    Opened files:\n"
       << "      fInput  = " << sInput.Data() << "\n"
       << "      fOutput = " << sOutput.Data()
       << endl;

  // grab input tuple
  TString sTupleToGrab = sTupleDir;
  sTupleToGrab.Append(sTuple.Data());

  TNtuple *ntToCalibrate = (TNtuple*) fInput -> Get(sTupleToGrab.Data());
  if (!ntToCalibrate) {
    cerr << "PANIC: couldn't grab input tuple!\n"
         << "       name  = " << sTupleToGrab << "\n"
         << "       tuple = " << ntToCalibrate << "\n"
         << endl;
    return;
  }
  cout << "    Grabbed input tuple:\n"
       << "      tuple = " << sTupleToGrab
       << endl;

  // declare tuple leaves
  Float_t ePar;
  Float_t fracParVsLeadBHCal;
  Float_t fracParVsLeadBEMC;
  Float_t fracParVsSumBHCal;
  Float_t fracParVsSumBEMC;
  Float_t fracLeadBHCalVsBEMC;
  Float_t fracSumBHCalVsBEMC;
  Float_t eLeadBHCal;
  Float_t eLeadBEMC;
  Float_t eSumBHCal;
  Float_t eSumBEMC;
  Float_t diffLeadBHCal;
  Float_t diffLeadBEMC;
  Float_t diffSumBHCal;
  Float_t diffSumBEMC;
  Float_t nHitsLeadBHCal;
  Float_t nHitsLeadBEMC;
  Float_t nClustBHCal;
  Float_t nClustBEMC;
  Float_t hLeadBHCal;
  Float_t hLeadBEMC;
  Float_t fLeadBHCal;
  Float_t fLeadBEMC;
  Float_t eLeadImage;
  Float_t eSumImage;
  Float_t eLeadSciFi;
  Float_t eSumSciFi;
  Float_t nClustImage;
  Float_t nClustSciFi;
  Float_t hLeadImage;
  Float_t hLeadSciFi;
  Float_t fLeadImage;
  Float_t fLeadSciFi;

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
  cout << "    Set tuple branches." << endl;

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
  const Float_t rEneBins[NRange]  = {-1.,   40.};
  const Float_t rDiffBins[NRange] = {-1.5,  5.5};
  const Float_t rFracBins[NRange] = {-0.05, 3.};

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
    hHCalEneBin[iEneBin]  = new TH1D(sHCalEne[iEneBin].Data(),  "", nEneBins,  rEneBins[0],  rEneBins[1]);
    hHCalDiffBin[iEneBin] = new TH1D(sHCalDiff[iEneBin].Data(), "", nDiffBins, rDiffBins[0], rDiffBins[1]);
    hHCalEneBin[iEneBin]  -> Sumw2();
    hHCalDiffBin[iEneBin] -> Sumw2();
  }
  cout << "    declared output histograms." << endl;

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
      const Bool_t isInEneParBin = ((ePar > eneParMin[iEneBin]) && (ePar < eneParMax[iEneBin]));
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
  Double_t  binSigmaEne[NEneBins];
  Double_t  valSigmaEne[NEneBins];
  Double_t  valSigmaDiff[NEneBins];
  Double_t  errSigmaEne[NEneBins];
  Double_t  errSigmaDiff[NEneBins];
  for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {

    // normalize hisotgrams
    const Double_t intEneBin  = hHCalEneBin[iEneBin]  -> Integral();
    const Double_t intDiffBin = hHCalDiffBin[iEneBin] -> Integral();
    if (intEneBin > 0.)  hHCalEneBin[iEneBin] -> Scale(1. / intEneBin);
    if (intDiffBin > 0.) hHCalDiffBin[iEneBin] -> Scale(1. / intDiffBin);

    // initialize functions
    fFitEneBin[iEneBin]  = new TF1(sFitEne[iEneBin].Data(),  "gaus(0)", xFitEneMin[iEneBin],  xFitEneMax[iEneBin]);
    fFitDiffBin[iEneBin] = new TF1(sFitDiff[iEneBin].Data(), "gaus(0)", xFitDiffMin[iEneBin], xFitDiffMax[iEneBin]);
    fFitEneBin[iEneBin]  -> SetParameter(0, ampEneGuess[iEneBin]);
    fFitEneBin[iEneBin]  -> SetParameter(1, muEneGuess[iEneBin]);
    fFitEneBin[iEneBin]  -> SetParameter(2, sigEneGuess[iEneBin]);
    fFitDiffBin[iEneBin] -> SetParameter(0, ampDiffGuess[iEneBin]);
    fFitDiffBin[iEneBin] -> SetParameter(1, muDiffGuess[iEneBin]);
    fFitDiffBin[iEneBin] -> SetParameter(2, sigDiffGuess[iEneBin]);
    fFitEneBin[iEneBin]  -> SetLineColor(fColEneBin[iEneBin]);
    fFitDiffBin[iEneBin] -> SetLineColor(fColEneBin[iEneBin]);

    // fit histograms
    hHCalEneBin[iEneBin]  -> Fit(sFitEne[iEneBin].Data(), "r");
    hHCalDiffBin[iEneBin] -> Fit(sFitDiff[iEneBin].Data(), "r");

    // grab resolutions and uncertainties
    const Double_t sigmaEne   = fFitEneBin[iEneBin]  -> GetParameter(2);
    const Double_t sigmaDiff  = fFitDiffBin[iEneBin] -> GetParameter(2);
    const Double_t errSigEne  = fFitEneBin[iEneBin]  -> GetParError(2);
    const Double_t errSigDiff = fFitDiffBin[iEneBin] -> GetParError(2);

    binSigmaEne[iEneBin]  = (eneParMin[iEneBin] - eneParMax[iEneBin]) / 2.;
    valSigmaEne[iEneBin]  = sigmaEne / enePar[iEneBin];
    valSigmaDiff[iEneBin] = sigmaDiff / enePar[iEneBin];
    errSigmaEne[iEneBin]  = errSigEne / enePar[iEneBin];
    errSigmaDiff[iEneBin] = errSigDiff / enePar[iEneBin];

    // set histogram styles
    hHCalEneBin[iEneBin]  -> SetMarkerColor(fColEneBin[iEneBin]);
    hHCalEneBin[iEneBin]  -> SetMarkerStyle(fMarEneBin[iEneBin]);
    hHCalEneBin[iEneBin]  -> SetLineColor(fColEneBin[iEneBin]);
    hHCalEneBin[iEneBin]  -> SetFillColor(fColEneBin[iEneBin]);
    hHCalEneBin[iEneBin]  -> GetXaxis() -> SetTitle(sEneTitleX.Data());
    hHCalDiffBin[iEneBin] -> SetMarkerColor(fColEneBin[iEneBin]);
    hHCalDiffBin[iEneBin] -> SetMarkerStyle(fMarEneBin[iEneBin]);
    hHCalDiffBin[iEneBin] -> SetLineColor(fColEneBin[iEneBin]);
    hHCalDiffBin[iEneBin] -> SetFillColor(fColEneBin[iEneBin]);
    hHCalDiffBin[iEneBin] -> GetXaxis() -> SetTitle(sDiffTitleX.Data());
  }
  cout << "    Normalized, fit, and set styles of resolution histograms." << endl;

  // create resolution graphs
  TGraphErrors *grResoEne  = new TGraphErrors(NEneBins, enePar, valSigmaEne,  binSigmaEne, errSigmaEne);
  TGraphErrors *grResoDiff = new TGraphErrors(NEneBins, enePar, valSigmaDiff, binSigmaEne, errSigmaDiff);
  grResoEne  -> SetName("grResoEne");
  grResoDiff -> SetName("grResoDiff");
  cout << "    Made resolution graphs." << endl;

  // plot fit distributions
  const UInt_t width(750);
  const UInt_t height(750);

  TCanvas *cResoEne = new TCanvas("cResoEne", "", width, height);
  cResoEne       -> cd();
  hHCalEneBin[0] -> Draw();
  for (UInt_t iEneBin = 1; iEneBin < NEneBins; iEneBin++) {
    hHCalEneBin[iEneBin] -> Draw("same");
  }
  fOutput  -> cd();
  cResoEne -> Write();
  cResoEne -> Close();

  TCanvas *cResoDiff = new TCanvas("cResoDiff", "", width, height);
  cResoDiff       -> cd();
  hHCalDiffBin[0] -> Draw();
  for (UInt_t iEneBin = 1; iEneBin < NEneBins; iEneBin++) {
    hHCalDiffBin[iEneBin] -> Draw("same");
  }
  fOutput   -> cd();
  cResoDiff -> Write();
  cResoDiff -> Close();
  cout << "    Made resolution plots." << endl;

  // do tmva training (if needed)
  TMVA::Factory    *factory;
  TMVA::DataLoader *loader;
  if (doTMVA) {

    // instantiate tmva library
    TMVA::Tools::Instance();
    cout << "    Beginning calibration:" << endl;

    // create tmva factory & load data
    factory = new TMVA::Factory("TMVARegression", fOutput, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");
    loader  = new TMVA::DataLoader("RegressionData");
    cout << "      Created factory and loaded data..." << endl;

    // set variables and target
    if (inSciGlassConfig) {
      if (addSpectators) {
        for (UInt_t iSpectator = 0; iSpectator < NSpecSci; iSpectator++) {
          loader -> AddSpectator(sSpecSci[iSpectator]);
        }
      }
      for (UInt_t iVariable = 0; iVariable < NVarSci; iVariable++) {
        loader -> AddVariable(sVarSci[iVariable]);
      }
    } else {
      if (addSpectators) {
        for (UInt_t iSpectator = 0; iSpectator < NSpecIma; iSpectator++) {
          loader -> AddSpectator(sSpecIma[iSpectator]);
        }
      }
      for (UInt_t iVariable = 0; iVariable < NVarIma; iVariable++) {
        loader -> AddVariable(sVarIma[iVariable]);
      }
    }
    loader -> AddTarget(sTarget);
    cout << "      Set spectators, variables, and target..." << endl;

    // add tree & prepare for training
    loader -> AddRegressionTree(ntToCalibrate, treeWeight);
    if (inSciGlassConfig) {
      loader -> PrepareTrainingAndTestTree(trainCutSci, "nTrain_Regression=1000:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V");
    } else {
      loader -> PrepareTrainingAndTestTree(trainCutIma, "nTrain_Regression=1000:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V");
    }
    cout << "      Added tree and prepared for training..." << endl;

    // book methods
    factory -> BookMethod(loader, TMVA::Types::kKNN, "KNN");
    factory -> BookMethod(loader, TMVA::Types::kLD,  "LD");
    factory -> BookMethod(loader, TMVA::Types::kMLP, "MLP");
    factory -> BookMethod(loader, TMVA::Types::kBDT, "BDTG");
    cout << "      Booked methods..." << endl;

    // train, test, & evaluate
    factory -> TrainAllMethods();
    factory -> TestAllMethods();
    factory -> EvaluateAllMethods();
    cout << "      Trained TMVA.\n"
         << "    Finished calibration!"
         << endl;
  }  // end if (doTmva)

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

  dReso      -> cd();
  grResoEne  -> Write();
  grResoDiff -> Write();
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
  if (doTMVA) {
    delete factory;
    delete loader;
  }
  return;

}

// end ------------------------------------------------------------------------
