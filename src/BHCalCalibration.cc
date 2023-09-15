// ----------------------------------------------------------------------------
// 'BHCalCalibration.cc'
// Derek Anderson
// 09.14.2023
//
// A simple interface to calibrate the simulated ePIC BHCal response with
// TMVA and calculate its energy resolution.
// ----------------------------------------------------------------------------

#define BHCALCALIBRATION_CC

// c utilities
#include <map>
#include <iostream>
// tmva classes
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"
// dataframe related classes
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>
// class header
#include "BHCalCalibration.hh"

using namespace std;
using namespace TMVA;



// public methods -------------------------------------------------------------

void BHCalCalibration::Init() {

  ParseInput();
  InitTuples();
  InitHistos();

  cout << "INIT" << endl;
  return;

}  // end 'Init()'



void BHCalCalibration::Train() {

  // instantiate tmva library
  Factory*    factory;
  DataLoader* loader;
  Tools::Instance();
  cout << "    Training TMVA:" << endl;

  // create tmva factory & load data
  factory = new Factory(_sFactory.data(), _fOutput, _sFactOpt.data());
  loader  = new DataLoader(_sLoader.data());
  cout << "      Created factory and loaded data." << endl;

  // add spectators if necessary
  if (_addSpecs) {
    for (const string spectator : _vecSpecForTMVA) {
      loader -> AddSpectator(spectator.data());
    }
    cout << "      Set spectators." << endl;
  }

  // add training variables
  for (const string variable : _vecVarsForTMVA) {
    loader -> AddVariable(variable.data());
  }
  cout << "      Set training variables." << endl;

  // set targets
  for (const string target : _vecTarsForTMVA) {
    loader -> AddTarget(target.data());
  }
  cout << "      Set regression targets." << endl;

  // add tree & prepare for training
  loader -> AddRegressionTree(_ntInput, _weight);
  loader -> PrepareTrainingAndTestTree(_cSelect, _sTrainOpt.data());
  cout << "      Added tree and prepared for training..." << endl;

  // book methods
  factory -> BookMethod(loader, Types::kLD,  "LD");
  factory -> BookMethod(loader, Types::kMLP, "MLP");
  factory -> BookMethod(loader, Types::kBDT, "BDTG");
  cout << "      Booked methods." << endl;

  // train, test, & evaluate
  factory -> TrainAllMethods();
  factory -> TestAllMethods();
  factory -> EvaluateAllMethods();
  cout << "      Trained TMVA!" << endl;
  return;

}  // end 'Train()'



void BHCalCalibration::Apply() {

  // default methods to be trained + tested
  map<string, int> Use;
  for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
    const string sToUse(sMethods[iMethod].data());
    Use[sToUse] = 1;
  }
  cout << "\n==> Start TMVARegressionApplication" << endl;

  Reader *reader = new Reader( "!Color:!Silent" );
  reader -> AddVariable("eLeadBHCal",       &eLeadBHCal);
  reader -> AddVariable("eLeadBEMC",        &eLeadBEMC);
  reader -> AddVariable("hLeadBHCal",       &hLeadBHCal);
  reader -> AddVariable("hLeadBEMC",        &hLeadBEMC);
  reader -> AddVariable("fLeadBHCal",       &fLeadBHCal);
  reader -> AddVariable("fLeadBEMC",        &fLeadBEMC);
  reader -> AddVariable("nHitsLeadBHCal",   &nHitsLeadBHCal);
  reader -> AddVariable("nHitsLeadBEMC",    &nHitsLeadBEMC);
  reader -> AddVariable("eSumImage",        &eSumImage);
  reader -> AddVariable("eSumSciFi",        &eSumSciFi);
  reader -> AddVariable("eSumSciFiLayer1",  &eSumSciFiLayer1);
  reader -> AddVariable("eSumSciFiLayer2",  &eSumSciFiLayer2);
  reader -> AddVariable("eSumSciFiLayer3",  &eSumSciFiLayer3);
  reader -> AddVariable("eSumSciFiLayer4",  &eSumSciFiLayer4);
  reader -> AddVariable("eSumSciFiLayer5",  &eSumSciFiLayer5);
  reader -> AddVariable("eSumSciFiLayer6",  &eSumSciFiLayer6);
  reader -> AddVariable("eSumSciFiLayer7",  &eSumSciFiLayer7);
  reader -> AddVariable("eSumSciFiLayer8",  &eSumSciFiLayer8);
  reader -> AddVariable("eSumSciFiLayer9",  &eSumSciFiLayer9);
  reader -> AddVariable("eSumSciFiLayer10", &eSumSciFiLayer10);
  reader -> AddVariable("eSumSciFiLayer11", &eSumSciFiLayer11);
  reader -> AddVariable("eSumSciFiLayer12", &eSumSciFiLayer12);
  reader -> AddVariable("eSumImageLayer1",  &eSumImageLayer1);
  reader -> AddVariable("eSumImageLayer2",  &eSumImageLayer2);
  reader -> AddVariable("eSumImageLayer3",  &eSumImageLayer3);
  reader -> AddVariable("eSumImageLayer4",  &eSumImageLayer4);
  reader -> AddVariable("eSumImageLayer5",  &eSumImageLayer5);
  reader -> AddVariable("eSumImageLayer6",  &eSumImageLayer6);

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

  // begin event loop
  TStopwatch stopwatch;
  cout << "--- Processing: " << nEvts << " events" << endl;

  nBytes = 0;
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
      const bool isInECalEneRange = ((eLeadBEMC > eneECalRange[0]) && (eLeadBEMC < eneECalRange[1]));
      if (doECalCut && !isInECalEneRange) continue;

      // fill resolution histograms
      if (methodExists) {
        for (UInt_t iCalibBin = 0; iCalibBin < NCalibBins; iCalibBin++) {
          const bool isInEneCalibBin = ((ePar > eneCalibMin[iCalibBin]) && (ePar < eneCalibMax[iCalibBin]));
          if (isInEneCalibBin) {
            hHCalCalibBin[method][iCalibBin] -> Fill(target);
          }
        }  // end energy bin loop
        hCalibCalibVsPar[method]  -> Fill(ePar,      target);
        hHCalCalibVsPar[method]   -> Fill(ePar,      eLeadBHCal);
        hHCalCalibVsCalib[method] -> Fill(target,    eLeadBHCal);
        hHCalCalibVsECal[method]  -> Fill(eLeadBEMC, eLeadBHCal);
        hECalCalibVsPar[method]   -> Fill(ePar,      eLeadBEMC);
        hECalCalibVsCalib[method] -> Fill(target,    eLeadBEMC);
      }  // end if (methodExists)
    }  // end method loop
  }  // end event loop
  stopwatch.Stop();

  // announce end of event loop
  cout << "--- End of event loop: ";
  cout << "\n    Application finished!" << endl;
  stopwatch.Print();
  return;

}  // end 'Apply()'



void BHCalCalibration::End() {

  FillHistos();
  ComputeReso();
  SaveOutput();

  cout << "END" << endl;
  return;

}  // end 'End()'



// private methods ------------------------------------------------------------

void BHCalCalibration::ParseInput() {

  cout << "  PARSE INPUT" << endl;
  return;

}  // end 'ParseInput()'



void BHCalCalibration::InitTuples() {

  cout << "  INIT TUPLES" << endl;
  return;

}  // end 'InitTuples()'



void BHCalCalibration::InitHistos() {

  cout << "  INIT HISTOS" << endl;
  return;

}  // end 'InitHistos()'



void BHCalCalibration::FillHistos() {

  cout << "  FILL HISTOS" << endl;
  return;

}  // end 'FillHistos()'



void BHCalCalibration::ComputeReso() {

  cout << "  COMPUTE RESO" << endl;
  return;

}  // end 'ComputeReso()'



void BHCalCalibration::SaveOutput() {

  cout << "  SAVE OUTPUT" << endl;
  return;

}  // end 'SaveOutput()'

// end ------------------------------------------------------------------------
