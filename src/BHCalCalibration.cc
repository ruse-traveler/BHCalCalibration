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
#include <cassert>
#include <iostream>
// dataframe related classes
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>
// class header
#include "BHCalCalibration.hh"

using namespace std;
using namespace TMVA;



// analysis methods -----------------------------------------------------------

void BHCalCalibration::Init() {

  OpenFiles();
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
    for (const string spectator : _vecSpecsForTMVA) {
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
  for (const string target : _vecTargsForTMVA) {
    loader -> AddTarget(target.data());
  }
  cout << "      Set regression targets." << endl;

  // add tree & prepare for training
  loader -> AddRegressionTree(_ntInput, _weight);
  loader -> PrepareTrainingAndTestTree(_cSelect, _sTrainOpt.data());
  cout << "      Added tree and prepared for training..." << endl;

  // book methods
  for (const auto methodAndOpts : _vecMethodsAndOptsTMVA) {
    if (get<1>(methodAndOpts).empty()) {
      factory -> BookMethod(loader, get<2>(methodAndOpts), get<0>(methodAndOpts).data());
    } else {
      factory -> BookMethod(loader, get<2>(methodAndOpts), get<0>(methodAndOpts).data(), get<1>(methodAndOpts).data());
    }
  }
  cout << "      Booked methods." << endl;

  // train, test, & evaluate
  factory -> TrainAllMethods();
  factory -> TestAllMethods();
  factory -> EvaluateAllMethods();
  cout << "      Trained TMVA!" << endl;

  delete factory;
  delete loader;
  return;

}  // end 'Train()'



void BHCalCalibration::Apply() {

  // default methods to be trained + tested
  map<string, int> Use;
  for (UInt_t iMethod = 0; iMethod < NMethods; iMethod++) {
    const string sToUse(sMethods[iMethod].data());
    Use[sToUse] = 1;
  }
  cout << "    Starting TMVA appllication:" << endl;

  // create reader
  Reader* reader = new Reader(_sReadOpt.data());
  for (const string variable : _vecVarsForTMVA) {
    reader -> AddVariable(variable.data(), &_mapInTupleVars[variable]);
  }
  cout << "      Created reader." << endl;

  // book methods to use
  for (const auto methodAndOpts : _vecMethodsAndOptsTMVA) {
    const string name    = get<0>(methodAndOpts);
    const string weights = _sLoader + "/weights/" + _sFactory + get<0>(methodAndOpts) + ".weights.xml";
    reader -> BookMVA(name, weights);
  }
  cout << "      Booked methods." << endl;

/* TODO  add TMVA histograms
  const UInt_t nTmvaBins(100);
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
*/

  // begin event loop
  const uint64_t nEntries = _ntInput -> GetEntries();
  cout << "      Beginning event loop:" << nEntries << " entries to process" << endl;

  uint64_t nBytes = 0;
  for (uint64_t iEntry = 0; iEntry < nEntries; iEntry++) {

    // grab entry
    const uint64_t bytes = _ntInput -> GetEntry(iEntry);
    if (bytes < 0.) {
      cerr << "WARNING: something wrong with entry " << iEntry << "! Aborting loop!" << endl;
      break;
    }
    nBytes += bytes;

    // announce progress
    const uint64_t iProg = iEntry + 1;
    if (iProg == nEntries) {
      cout << "        Processing entry " << iProg << "/" << nEntries << "..." << endl;
    } else {
      cout << "        Processing entry " << iProg << "/" << nEntries << "...\r" << flush;
    }

    // loop over methods
    for (const auto methodAndOpt : _vecMethodsAndOptsTMVA) {

      // grab regression target
      const float target = (reader -> EvaluateRegression(get<0>(methodAndOpt)))[0];

      // TODO fill TMVA hist
      //TString title  = hTMVA[iTmvaHist] -> GetTitle();
      //hTMVA[iTmvaHist] -> Fill(target);

      /* TODO fill output histograms */

    }  // end method loop
  }  // end event loop

  // announce end of event loop
  cout << "      Application finished!" << endl;

  // delete reader and exit
  delete reader;
  return;

}  // end 'Apply()'



void BHCalCalibration::End() {

  FillHistos();
  ComputeReso();
  SaveOutput();

  cout << "END" << endl;
  return;

}  // end 'End()'



// setters --------------------------------------------------------------------

void BHCalCalibration::SetInput(const string sInFile, const string sInTuple, const float treeWeight) {

  _sInFile  = sInFile;
  _sInTuple = sInTuple;
  _weight   = treeWeight;

  cout << "SET INPUT" << endl;
  return;

}  // end 'SetInput(string, string, float)'



void BHCalCalibration::SetTupleArgs(const vector<string> vecInput) {

  _vecInTupleLeaves = vecInput;
  for (const string inputLeaf : _vecInTupleLeaves) {
    _mapInTupleVars[inputLeaf] = -999.;
  }

  cout << "SET TUPLE ARGS" << endl;
  return;

}  // end 'SetTupleArgs(vector<string>)'



void BHCalCalibration::SetTmvaOpts(const string sFactOpt, const string sTrainOpt, const string sReadOpt, const bool addSpecs) {

  _sFactOpt  = sFactOpt;
  _sTrainOpt = sTrainOpt;
  _sReadOpt  = sReadOpt;
  _addSpecs  = addSpecs;

  cout << "SET TMVA OPTS" << endl;
  return;

}  // end 'SetTmvaOpts(string, string, string, bool)'



void BHCalCalibration::SetTmvaArgs(const vector<string> vecVars, const vector<string> vecTargs, const vector<string> vecSpecs = {}, const TCut select = "") {

  _vecVarsForTMVA  = vecVars;
  _vecTargsForTMVA = vecTargs;
  _vecSpecsForTMVA = vecSpecs;
  _cSelect         = select;

  cout << "SET TMVA ARGS" << endl;
  return;

}  // end 'SetTmvaArgs(vector<string>, vector<string>, vector<string>, TCut)'



void BHCalCalibration::SetTmvaMethods(const vector<pair<string, string>> vecMethodAndOpts) {

  for (const auto methodAndOpt : vecMethodAndOpts) {
    _vecMethodsAndOptsTMVA.push_back(make_tuple(methodAndOpt.first, methodAndOpt.second, _mapMethodToIndex[methodAndOpt.first]));
  }

  /* TODO check for duplicate methods */

  cout << "SET TMVA METHODS" << endl;
  return;

}  // end 'SetTmvaMethods(pair<vector<string, vector<string>>)'



// private methods ------------------------------------------------------------

void BHCalCalibration::OpenFiles() {

  _fInput  = new TFile(_sInFile.data(),  "read");
  _fOutput = new TFile(_sOutFile.data(), "recreate");
  if (!_fInput || !_fOutput) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fInput = " << _fInput << ", fOutput = " << _fOutput << "\n"
         << endl;
    assert(_fInput && _fOutput);
  }

  cout << "  OPEN FILES" << endl;
  return;

}  // end 'OpenFiles()'


void BHCalCalibration::InitTuples() {

  // grab input tuple
  _ntInput = (TNtuple*) _fInput -> Get(_sInTuple.data());
  if (!_ntInput) {
    cerr << "PANIC: couldn't grab input tuple!\n"
         << "       tuple = " << _sInTuple << " = " << _ntInput << "\n"
         << endl;
    assert(_ntInput);
  }

  // set input branches
  for (const string leaf : _vecInTupleLeaves) {
    _ntInput -> SetBranchAddress(leaf.data(), &_mapInTupleVars[leaf]);
  }

  /* TODO set up output tuple */

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



void BHCalCalibration::CloseFiles() {

  _fOutput -> cd();
  _fOutput -> Close();
  _fInput  -> cd();
  _fInput  -> Close();

  cout << "  CLOSE FILES" << endl;
  return;

}  // end 'CloseFiles()'

// end ------------------------------------------------------------------------
