// ----------------------------------------------------------------------------
// 'ApplyHCalCalibration.C'
// Derek Anderson
// 03.28.2023
//
// Use this to apply the TMVA training done
// in 'DoHCalCalibration.C'.
//
// Note: 'fConfig' sets which BEMC configuration
// is being used.
//   fConfig = 0: SciGlass BEMC
//   fConfig = 1: Imaging BEMC
//   fConfig = 2: Default (see "parse configuration" block)
//
// Derived from code by Andreas Hoecker
// ---------------------------------------------------------------------------

// standard c includes
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
// root includes
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TDirectory.h"
// tmva includes
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace std;
using namespace TMVA;

// global constants
static const UInt_t NTxt(2);
static const UInt_t NVtx(4);
static const UInt_t NRange(2);
static const UInt_t NEneBins(4);
static const UInt_t NMethods(3);

// tmva constants
static const UInt_t  NTmvaHistMax(100);
static const TString STmvaPrefix("TMVARegression");
static const TString STmvaDirSci("tmva/SciGlassRegressionData_NoNClustAndWithNHits/weights/");
static const TString STmvaDirIma("tmva/ImagingRegressionData_NoNClustAndWithNHits/weights/");

// default arguments
static const UInt_t  FConfigDef(0);
static const TString SInSciDef("./eicrecon_output/merged/forECalStudy.sciglass.e2t20th35145n5KeaPim.d8m3y2023.plugin.root");
static const TString SInImaDef("./eicrecon_output/merged/forECalStudy.imaging.e2t20th35145n5KeaPim.d8m3y2023.plugin.root");
static const TString SOutSciDef("forSciGlassReso.application_forMipCheck_ecalEneG05.e2t20th35145n5KeaPim.d14m3y2023.root");
static const TString SOutImaDef("forImagingReso.application_forMipCheck_ecalEneG05.e2t20th35145n5KeaPim.d14m3y2023.root");
static const TString STupleDef("ntForCalibration");



void ApplyHCalCalibration(const UInt_t fConfig = FConfigDef, const TString sInSci = SInSciDef, const TString sInIma = SInImaDef, const TString sOutSci = SOutSciDef, const TString sOutIma = SOutImaDef, const TString sTuple = STupleDef, const string myMethodList = "") {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning BHCal calibration..." << endl; 

  // ecal cut parameters
  const Bool_t   doECalCut(true);
  const Double_t eneECalRange[NRange] = {0.5, 100.};

  // histogram parameters
  const UInt_t  fColEneBin[NEneBins]   = {809, 909, 889, 869};
  const UInt_t  fMarEneBin[NEneBins]   = {26,  27,  24,  25};
  const TString sHCalEneBase[NEneBins] = {"hHCalEne_ene2",  "hHCalEne_ene5",  "hHCalEne_ene10",  "hHCalEne_ene20"};
  const TString sMethods[NMethods]     = {"LD", "MLP", "BDTG"};
  const TString sEneTitleX("E_{par}^{reco} [GeV]");
  const TString sTitleY("arbitrary units");

  // generic resolution parameters
  const Double_t enePar[NEneBins]       = {2., 5., 10., 20.};
  const Double_t eneParMin[NEneBins]    = {1., 3., 7.,  13.};
  const Double_t eneParMax[NEneBins]    = {3., 7., 13., 27.};

  // reco vs. par ene resolution parameters
  const Double_t xFitEneMin[NEneBins]   = {0.5, 4.,  8.,  13.};
  const Double_t xFitEneMax[NEneBins]   = {5.5, 8.,  14., 23.};
  const Double_t ampEneGuess[NEneBins]  = {1.,  1.,  1.,  1.};
  const Double_t muEneGuess[NEneBins]   = {3.,  6.,  11., 18.};
  const Double_t sigEneGuess[NEneBins]  = {2.,  2.,  3.,  5.};
  const TString  sFitEneBase[NEneBins]  = {"fFitEne_ene2", "fFitEne_ene5", "fFitEne_ene10", "fFitEne_ene20"};

  // style parameters
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCenter(1);
  const Float_t fOffX(1.2);
  const Float_t fOffY(1.3);
  const TString sTitle("");

  // text parameters
  const TString sHeader("");
  const TString sTxt[NTxt]       = {"ePIC simulation [23.01.0]", "single #pi^{-}"};
  const TString sLabel[NEneBins] = {"E_{par} = 2 GeV", "E_{par} = 5 GeV", "E_{par} = 10 GeV", "E_{par} = 20 GeV"};

  // parse configuration
  Bool_t  isInSciGlassConfig;
  TString sInput;
  TString sOutput;
  TString sTupleDir;
  TString sTmvaDir;
  switch (fConfig) {
    case 0:
      isInSciGlassConfig = true;
      sInput             = sInSci;
      sOutput            = sOutSci;
      sTupleDir          = "JCalibrateHCalWithSciGlass/";
      sTmvaDir           = STmvaDirSci;
      break;
    case 1:
      isInSciGlassConfig = false;
      sInput             = sInIma;
      sOutput            = sOutIma;
      sTupleDir          = "JCalibrateHCalWithImaging/";
      sTmvaDir           = STmvaDirIma;
      break;
    default:
      isInSciGlassConfig = true;
      sInput             = sInSci;
      sOutput            = sOutSci;
      sTupleDir          = "JCalibrateHCal/";
      sTmvaDir           = STmvaDirSci;
      break;
  }

  if (isInSciGlassConfig) {
    cout << "    Using SciGlass configuration..." << endl;
  } else {
    cout << "    Using imaging configuration..." << endl;
  }

  // open files
  TFile *fInput  = new TFile(sInput.Data(),  "read");
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

  // tuple leaves
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

  // resolution histograms
  TH1D *hHCalEneBin[NMethods][NEneBins];
  TH2D *hCalibEneVsPar[NMethods];
  TH2D *hHCalEneVsPar[NMethods];
  TH2D *hHCalEneVsCalib[NMethods];
  TH2D *hHCalEneVsECal[NMethods];
  TH2D *hECalEneVsPar[NMethods];
  TH2D *hECalEneVsCalib[NMethods];

  // histogram binning
  const Ssiz_t  nEneBins(41);
  const Ssiz_t  nEneBins2D(410);
  const Ssiz_t  nDiffBins(700);
  const Ssiz_t  nFracBins(305);
  const Float_t rEneBins[NRange]  = {-1.,   40.};
  const Float_t rDiffBins[NRange] = {-1.5,  5.5};
  const Float_t rFracBins[NRange] = {-0.05, 3.};

  // declare resolution histograms
  for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
    for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {

      // make name
      TString sHCalEne(sHCalEneBase[iEneBin].Data());
      sHCalEne.Append("_");
      sHCalEne.Append(sMethods[iMethod].Data());
     
      hHCalEneBin[iMethod][iEneBin] = new TH1D(sHCalEne.Data(), "", nEneBins, rEneBins[0], rEneBins[1]);
      hHCalEneBin[iMethod][iEneBin] -> Sumw2();
    }

    // make name
    TString sCalibEneVsPar("hCalibEneVsPar");
    TString sHCalEneVsPar("hHCalEneVsPar");
    TString sHCalEneVsCalib("hHCalEneVsCalib");
    TString sHCalEneVsECal("hHCalEneVsECal");
    TString sECalEneVsPar("hECalEneVsPar");
    TString sECalEneVsCalib("hECalEneVsCalib");
    sCalibEneVsPar.Append("_");
    sHCalEneVsPar.Append("_");
    sHCalEneVsCalib.Append("_");
    sHCalEneVsECal.Append("_");
    sECalEneVsPar.Append("_");
    sECalEneVsCalib.Append("_");
    sCalibEneVsPar.Append(sMethods[iMethod].Data());
    sHCalEneVsPar.Append(sMethods[iMethod].Data());
    sHCalEneVsCalib.Append(sMethods[iMethod].Data());
    sHCalEneVsECal.Append(sMethods[iMethod].Data());
    sECalEneVsPar.Append(sMethods[iMethod].Data());
    sECalEneVsCalib.Append(sMethods[iMethod].Data());

    hCalibEneVsPar[iMethod]  = new TH2D(sCalibEneVsPar.Data(),  "", nEneBins,   rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hHCalEneVsPar[iMethod]   = new TH2D(sHCalEneVsPar.Data(),   "", nEneBins2D, rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hHCalEneVsCalib[iMethod] = new TH2D(sHCalEneVsCalib.Data(), "", nEneBins2D, rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hHCalEneVsECal[iMethod]  = new TH2D(sHCalEneVsECal.Data(),  "", nEneBins2D, rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hECalEneVsPar[iMethod]   = new TH2D(sECalEneVsPar.Data(),   "", nEneBins2D, rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hECalEneVsCalib[iMethod] = new TH2D(sECalEneVsCalib.Data(), "", nEneBins2D, rEneBins[0], rEneBins[1], nEneBins, rEneBins[0], rEneBins[1]);
    hCalibEneVsPar[iMethod]  -> Sumw2();
    hHCalEneVsPar[iMethod]   -> Sumw2();
    hHCalEneVsCalib[iMethod] -> Sumw2();
    hHCalEneVsECal[iMethod]  -> Sumw2();
    hECalEneVsPar[iMethod]   -> Sumw2();
    hECalEneVsCalib[iMethod] -> Sumw2();
  }  // end method loop
  cout << "    Declared resolution histograms.\n"
       << "    Beginning application..."
       << endl;

  // instantiate tmva library
  Tools::Instance();

  // default methods to be trained + tested
  map<string, int> Use;
  for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
    const string sToUse(sMethods[iMethod].Data());
    Use[sToUse] = 1;
  }
  cout << "\n==> Start TMVARegressionApplication" << endl;

  Reader *reader = new Reader( "!Color:!Silent" );
  reader -> AddVariable("eLeadBHCal",     &eLeadBHCal);
  reader -> AddVariable("eLeadBEMC",      &eLeadBEMC);
  reader -> AddVariable("hLeadBHCal",     &hLeadBHCal);
  reader -> AddVariable("hLeadBEMC",      &hLeadBEMC);
  reader -> AddVariable("fLeadBHCal",     &fLeadBHCal);
  reader -> AddVariable("fLeadBEMC",      &fLeadBEMC);
  reader -> AddVariable("nHitsLeadBHCal", &nHitsLeadBHCal);
  reader -> AddVariable("nHitsLeadBEMC",  &nHitsLeadBEMC);
  if (!isInSciGlassConfig) {
    reader -> AddVariable("eSumImage",    &eSumImage);
    reader -> AddVariable("eSumSciFi",    &eSumSciFi);
  }

  // book method(s)
  for (map<string, int>::iterator itMethod = Use.begin(); itMethod != Use.end(); itMethod++) {
    if (itMethod -> second) {
      TString methodName = TString(itMethod -> first) + " method";
      TString weightfile = sTmvaDir + STmvaPrefix + "_" + TString(itMethod -> first) + ".weights.xml";
      reader->BookMVA(methodName, weightfile);
    }
  }  // end method loop

  // for tmva histogram binning
  const UInt_t  nTmvaBins(100);
  const Float_t rTmvaBins[NRange] = {-100., 600.};

  // Book tmva histograms
  Int_t  nTmvaHist(-1);
  TH1   *hTMVA[NTmvaHistMax];
  for (map<string, int>::iterator itMethod = Use.begin(); itMethod != Use.end(); itMethod++) {
    TString  sName  = TString(itMethod -> first.c_str());
    TString  sTitle = TString(itMethod -> first) + " method";
    TH1     *hNew   = new TH1F(sName.Data(), sTitle.Data(), nTmvaBins, rTmvaBins[0], rTmvaBins[1]);
    if (!hNew) {
      cerr << "PANIC: couldn't create TMVA histogram #" << nTmvaHist << "! Aborting code execution!\n" << endl;
      return;
    } else {
      if (itMethod -> second) hTMVA[++nTmvaHist] = hNew;
    }
  }  // end method loop
  nTmvaHist++;

  // begin event loop
  TStopwatch stopwatch;
  Long64_t   nBytes = 0;
  Long64_t   nEvts  = ntToCalibrate -> GetEntries();
  cout << "--- Processing: " << nEvts << " events" << endl;

  stopwatch.Start();
  for (Long64_t iEvt = 0; iEvt < nEvts; iEvt++) {

    // announce progress
    if (iEvt % 1000 == 0) {
      cout << "--- ... Processing event: " << iEvt << endl;
    }

    const Long64_t bytes = ntToCalibrate -> GetEntry(iEvt);
    if (bytes < 0.) {
      cerr << "WARNING something wrong with event " << iEvt << "! Aborting loop!" << endl;
      break;
    }
    nBytes += bytes;

    // loop over methods
    for (Int_t iTmvaHist = 0; iTmvaHist < nTmvaHist; iTmvaHist++) {

      // grab regression target
      TString title  = hTMVA[iTmvaHist] -> GetTitle();
      Float_t target = (reader -> EvaluateRegression(title))[0];
      hTMVA[iTmvaHist] -> Fill(target);

      // check for method
      Int_t method = -1;
      for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
        Bool_t isMethod = title.Contains(sMethods[iMethod].Data());
        if (isMethod) {
          method = iMethod;
          break;
        }
      }  // end method loop

      // check for ecal energy
      const Bool_t methodExists     = (method > -1);
      const Bool_t isInECalEneRange = ((eLeadBEMC > eneECalRange[0]) && (eLeadBEMC < eneECalRange[1]));
      if (doECalCut && !isInECalEneRange) continue;

      // fill resolution histograms
      if (methodExists) {
        for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {
          const Bool_t isInEneParBin = ((ePar > eneParMin[iEneBin]) && (ePar < eneParMax[iEneBin]));
          if (isInEneParBin) {
            hHCalEneBin[method][iEneBin] -> Fill(target);
          }
        }  // end energy bin loop
        hCalibEneVsPar[method]  -> Fill(ePar,      target);
        hHCalEneVsPar[method]   -> Fill(ePar,      eLeadBHCal);
        hHCalEneVsCalib[method] -> Fill(target,    eLeadBHCal);
        hHCalEneVsECal[method]  -> Fill(eLeadBEMC, eLeadBHCal);
        hECalEneVsPar[method]   -> Fill(ePar,      eLeadBEMC);
        hECalEneVsCalib[method] -> Fill(target,    eLeadBEMC);
      }  // end if (methodExists)
    }  // end method loop
  }  // end event loop
  stopwatch.Stop();

  // announce end of event loop
  cout << "--- End of event loop: ";
  stopwatch.Print();

  // for graphs
  Double_t binSigmaEne[NMethods][NEneBins];
  Double_t valMuEne[NMethods][NEneBins];
  Double_t valMuEneHist[NMethods][NEneBins];
  Double_t valSigmaEne[NMethods][NEneBins];
  Double_t valSigmaEneHist[NMethods][NEneBins];
  Double_t errMuEne[NMethods][NEneBins];
  Double_t errMuEneHist[NMethods][NEneBins];
  Double_t errSigmaEne[NMethods][NEneBins];
  Double_t errSigmaEneHist[NMethods][NEneBins];
  cout << "\n    Application finished!" << endl;

  // resolution calculation
  TF1          *fFitEneBin[NMethods][NEneBins];
  TCanvas      *cResoEne[NMethods];
  TGraphErrors *grLineEne[NMethods];
  TGraphErrors *grLineEneHist[NMethods];
  TGraphErrors *grResoEne[NMethods];
  TGraphErrors *grResoEneHist[NMethods];
  for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
    for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {

      // normalize hisotgrams
      const Double_t intEneBin = hHCalEneBin[iMethod][iEneBin]  -> Integral();
      if (intEneBin > 0.) hHCalEneBin[iMethod][iEneBin] -> Scale(1. / intEneBin);

      // make name
      TString sFitEne(sFitEneBase[iEneBin].Data());
      sFitEne.Append("_");
      sFitEne.Append(sMethods[iMethod]);

      // initialize functions
      fFitEneBin[iMethod][iEneBin] = new TF1(sFitEne.Data(),  "gaus(0)", xFitEneMin[iEneBin],  xFitEneMax[iEneBin]);
      fFitEneBin[iMethod][iEneBin] -> SetParameter(0, ampEneGuess[iEneBin]);
      fFitEneBin[iMethod][iEneBin] -> SetParameter(1, muEneGuess[iEneBin]);
      fFitEneBin[iMethod][iEneBin] -> SetParameter(2, sigEneGuess[iEneBin]);
      fFitEneBin[iMethod][iEneBin] -> SetLineColor(fColEneBin[iEneBin]);

      // fit histograms
      hHCalEneBin[iMethod][iEneBin] -> Fit(sFitEne.Data(), "r");

      // grab resolutions and uncertainties
      const Double_t mu     = fFitEneBin[iMethod][iEneBin]  -> GetParameter(1);
      const Double_t sigma  = fFitEneBin[iMethod][iEneBin]  -> GetParameter(2);
      const Double_t errMu  = fFitEneBin[iMethod][iEneBin]  -> GetParError(1);
      const Double_t errSig = fFitEneBin[iMethod][iEneBin]  -> GetParError(2);
      const Double_t perMu  = errMu / mu;
      const Double_t perSig = errSig / sigma;

      const Double_t muHist     = hHCalEneBin[iMethod][iEneBin]  -> GetMean();
      const Double_t sigmaHist  = hHCalEneBin[iMethod][iEneBin]  -> GetRMS();
      const Double_t errMuHist  = hHCalEneBin[iMethod][iEneBin]  -> GetMeanError();
      const Double_t errSigHist = hHCalEneBin[iMethod][iEneBin]  -> GetRMSError();
      const Double_t perMuHist  = errMuHist / muHist;
      const Double_t perSigHist = errSigHist / sigmaHist;

      // set fit values
      binSigmaEne[iMethod][iEneBin] = (eneParMin[iEneBin] - eneParMax[iEneBin]) / 2.;
      valMuEne[iMethod][iEneBin]    = mu;
      valSigmaEne[iMethod][iEneBin] = sigma / mu;
      errMuEne[iMethod][iEneBin]    = errMu;
      errSigmaEne[iMethod][iEneBin] = valSigmaEne[iMethod][iEneBin] * TMath::Sqrt((perMu * perMu) + (perSig * perSig));

      // set histogram values
      valMuEneHist[iMethod][iEneBin]    = muHist;
      valSigmaEneHist[iMethod][iEneBin] = sigmaHist / muHist;
      errMuEneHist[iMethod][iEneBin]    = errMuHist;
      errSigmaEneHist[iMethod][iEneBin] = valSigmaEneHist[iMethod][iEneBin] * TMath::Sqrt((perMuHist * perMuHist) + (perSigHist * perSigHist));

      // set histogram styles
      hHCalEneBin[iMethod][iEneBin] -> SetMarkerColor(fColEneBin[iEneBin]);
      hHCalEneBin[iMethod][iEneBin] -> SetMarkerStyle(fMarEneBin[iEneBin]);
      hHCalEneBin[iMethod][iEneBin] -> SetLineColor(fColEneBin[iEneBin]);
      hHCalEneBin[iMethod][iEneBin] -> SetLineStyle(fLin);
      hHCalEneBin[iMethod][iEneBin] -> SetFillColor(fColEneBin[iEneBin]);
      hHCalEneBin[iMethod][iEneBin] -> SetFillStyle(fFil);
      hHCalEneBin[iMethod][iEneBin] -> SetTitle(sTitle.Data());
      hHCalEneBin[iMethod][iEneBin] -> SetTitleFont(fTxt);
      hHCalEneBin[iMethod][iEneBin] -> GetXaxis() -> SetTitle(sEneTitleX.Data());
      hHCalEneBin[iMethod][iEneBin] -> GetXaxis() -> SetTitleFont(fTxt);
      hHCalEneBin[iMethod][iEneBin] -> GetXaxis() -> SetTitleOffset(fOffX);
      hHCalEneBin[iMethod][iEneBin] -> GetXaxis() -> CenterTitle(fCenter);
      hHCalEneBin[iMethod][iEneBin] -> GetYaxis() -> SetTitle(sTitleY.Data());
      hHCalEneBin[iMethod][iEneBin] -> GetYaxis() -> SetTitleFont(fTxt);
      hHCalEneBin[iMethod][iEneBin] -> GetYaxis() -> SetTitleOffset(fOffY);
      hHCalEneBin[iMethod][iEneBin] -> GetYaxis() -> CenterTitle(fCenter);
    }
    cout << "    Fit resolution histograms and set styles." << endl;

    // make name
    TString sGraphLineEne("grLineEne");
    TString sGraphLineEneHist("grLineEneHist");
    TString sGraphResoEne("grResoEne");
    TString sGraphResoEneHist("grResoEneHist");
    sGraphLineEne.Append("_");
    sGraphLineEneHist.Append("_");
    sGraphResoEne.Append("_");
    sGraphResoEneHist.Append("_");
    sGraphLineEne.Append(sMethods[iMethod].Data());
    sGraphLineEneHist.Append(sMethods[iMethod].Data());
    sGraphResoEne.Append(sMethods[iMethod].Data());
    sGraphResoEneHist.Append(sMethods[iMethod].Data());

    // create resolution graphs
    grLineEne[iMethod]     = new TGraphErrors(NEneBins, enePar, valMuEne[iMethod],        binSigmaEne[iMethod], errMuEne[iMethod]);
    grLineEneHist[iMethod] = new TGraphErrors(NEneBins, enePar, valMuEneHist[iMethod],    binSigmaEne[iMethod], errMuEneHist[iMethod]);
    grResoEne[iMethod]     = new TGraphErrors(NEneBins, enePar, valSigmaEne[iMethod],     binSigmaEne[iMethod], errSigmaEne[iMethod]);
    grResoEneHist[iMethod] = new TGraphErrors(NEneBins, enePar, valSigmaEneHist[iMethod], binSigmaEne[iMethod], errSigmaEneHist[iMethod]);
    grLineEne[iMethod]     -> SetName(sGraphLineEne.Data());
    grLineEneHist[iMethod] -> SetName(sGraphLineEneHist.Data());
    grResoEne[iMethod]     -> SetName(sGraphResoEne.Data());
    grResoEneHist[iMethod] -> SetName(sGraphResoEneHist.Data());

    // make legend
    const UInt_t  fColLeg      = 0;
    const UInt_t  fFilLeg      = 0;
    const UInt_t  fLinLeg      = 0;
    const Float_t hObjLeg      = NEneBins * 0.05;
    const Float_t yObjLeg      = 0.1 + hObjLeg;
    const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, yObjLeg};

    TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sHeader.Data());
    leg -> SetFillColor(fColLeg);
    leg -> SetFillStyle(fFilLeg);
    leg -> SetLineColor(fColLeg);
    leg -> SetLineStyle(fLinLeg);
    leg -> SetTextFont(fTxt);
    leg -> SetTextAlign(fAln);
    for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {
      leg -> AddEntry(hHCalEneBin[iMethod][iEneBin], sLabel[iEneBin].Data(), "pf");
    }
    cout << "    Made legend." << endl;

    // make text
    const UInt_t  fColTxt      = 0;
    const UInt_t  fFilTxt      = 0;
    const UInt_t  fLinTxt      = 0;
    const Float_t hObjTxt      = NTxt * 0.05;
    const Float_t yObjTxt      = 0.1 + hObjTxt;
    const Float_t fTxtXY[NVtx] = {0.3, 0.1, 0.5, yObjTxt};

    TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
    txt -> SetFillColor(fColTxt);
    txt -> SetFillStyle(fFilTxt);
    txt -> SetLineColor(fColTxt);
    txt -> SetLineStyle(fLinTxt);
    txt -> SetTextFont(fTxt);
    txt -> SetTextAlign(fAln);
    for (UInt_t iTxt = 0; iTxt < NTxt; iTxt++) {
      txt -> AddText(sTxt[iTxt].Data());
    }
    cout << "    Made text." << endl;

    // plot fit distributions
    const UInt_t  width(750);
    const UInt_t  height(750);
    const UInt_t  fMode(0);
    const UInt_t  fBord(2);
    const UInt_t  fGrid(0);
    const UInt_t  fTick(1);
    const UInt_t  fLogX(0);
    const UInt_t  fLogY(1);
    const UInt_t  fFrame(0);
    const Float_t fMarginL(0.15);
    const Float_t fMarginR(0.02);
    const Float_t fMarginT(0.02);
    const Float_t fMarginB(0.15);

    // make name
    TString sResoEne("cResoEne");
    sResoEne.Append("_");
    sResoEne.Append(sMethods[iMethod].Data());

    cResoEne[iMethod] = new TCanvas(sResoEne.Data(), "", width, height);
    cResoEne[iMethod]       -> SetGrid(fGrid, fGrid);
    cResoEne[iMethod]       -> SetTicks(fTick, fTick);
    cResoEne[iMethod]       -> SetBorderMode(fMode);
    cResoEne[iMethod]       -> SetBorderSize(fBord);
    cResoEne[iMethod]       -> SetFrameBorderMode(fFrame);
    cResoEne[iMethod]       -> SetLeftMargin(fMarginL);
    cResoEne[iMethod]       -> SetRightMargin(fMarginR);
    cResoEne[iMethod]       -> SetTopMargin(fMarginT);
    cResoEne[iMethod]       -> SetBottomMargin(fMarginB);
    cResoEne[iMethod]       -> SetLogx(fLogX);
    cResoEne[iMethod]       -> SetLogy(fLogY);
    cResoEne[iMethod]       -> cd();
    hHCalEneBin[iMethod][0] -> Draw();
    for (UInt_t iEneBin = 1; iEneBin < NEneBins; iEneBin++) {
      hHCalEneBin[iMethod][iEneBin] -> Draw("same");
    }
    leg               -> Draw();
    txt               -> Draw();
    fOutput           -> cd();
    cResoEne[iMethod] -> Write();
    cResoEne[iMethod] -> Close();
  }  // end method loop
  cout << "    Made resolution plots." << endl;

  // create directories
  TDirectory *dReso = (TDirectory*) fOutput -> mkdir("reso");
  TDirectory *dTmva = (TDirectory*) fOutput -> mkdir("tmva");

  // write histograms
  dReso -> cd();
  for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
    hCalibEneVsPar[iMethod]  -> Write();
    hHCalEneVsPar[iMethod]   -> Write();
    hHCalEneVsCalib[iMethod] -> Write();
    hHCalEneVsECal[iMethod]  -> Write();
    hECalEneVsPar[iMethod]   -> Write();
    hECalEneVsCalib[iMethod] -> Write();
    grLineEne[iMethod]       -> Write();
    grLineEneHist[iMethod]   -> Write();
    grResoEne[iMethod]       -> Write();
    grResoEneHist[iMethod]   -> Write();
    for (UInt_t iEneBin = 0; iEneBin < NEneBins; iEneBin++) {
      hHCalEneBin[iMethod][iEneBin] -> Write();
      fFitEneBin[iMethod][iEneBin]  -> Write();
    }
  }  // end method loop

  // write tmva histograms
  dTmva -> cd();
  for (Int_t iTmvaHist = 0; iTmvaHist < nTmvaHist; iTmvaHist++) {
    hTMVA[iTmvaHist] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();

  // delete tmva reader and exit
  delete reader;
  cout << "  Finished Calibration application script!\n" << endl;

}

// end ------------------------------------------------------------------------
