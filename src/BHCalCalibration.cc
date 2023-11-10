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
#include <TDirectory.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>
// class header
#include "BHCalCalibration.hh"

using namespace std;
using namespace TMVA;



// analysis methods -----------------------------------------------------------

void BHCalCalibration::Init() {

  // announce initialization
  cout << "\n  Beginning BHCal calibration...\n"
       << "    Initializing:"
       << endl;

  OpenFiles();
  InitTuples();
  InitHistos();
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
    switch (get<2>(methodAndOpts)) {
      case TMVA::Types::kPDERS:
        factory -> BookMethod(loader, TMVA::Types::kPDERS, get<0>(methodAndOpts).data(), get<1>(methodAndOpts).data());
        break;
      case TMVA::Types::kPDEFoam:
        factory -> BookMethod(loader, TMVA::Types::kPDEFoam, get<0>(methodAndOpts).data(), get<1>(methodAndOpts).data());
        break;
      case TMVA::Types::kKNN:
        factory -> BookMethod(loader, TMVA::Types::kKNN, get<0>(methodAndOpts).data(), get<1>(methodAndOpts).data());
        break;
      case TMVA::Types::kLD:
        factory -> BookMethod(loader, TMVA::Types::kLD, get<0>(methodAndOpts).data(), get<1>(methodAndOpts).data());
        break;
      case TMVA::Types::kFDA:
        factory -> BookMethod(loader, TMVA::Types::kFDA, get<0>(methodAndOpts).data(), get<1>(methodAndOpts).data());
        break;
      case TMVA::Types::kMLP:
        factory -> BookMethod(loader, TMVA::Types::kMLP, get<0>(methodAndOpts).data(), get<1>(methodAndOpts).data());
        break;
      case TMVA::Types::kDL:
        factory -> BookMethod(loader, TMVA::Types::kDL, get<0>(methodAndOpts).data(), get<1>(methodAndOpts).data());
        break;
      case TMVA::Types::kSVM:
        factory -> BookMethod(loader, TMVA::Types::kSVM, get<0>(methodAndOpts).data(), get<1>(methodAndOpts).data());
        break;
      case TMVA::Types::kBDT:
        factory -> BookMethod(loader, TMVA::Types::kBDT, get<0>(methodAndOpts).data(), get<1>(methodAndOpts).data());
        break;
    }
  }
  cout << "      Booked methods." << endl;

  // train, test, & evaluate
  factory -> TrainAllMethods();
  factory -> TestAllMethods();
  factory -> EvaluateAllMethods();
  cout << "    Trained TMVA!" << endl;

  delete factory;
  delete loader;
  return;

}  // end 'Train()'



void BHCalCalibration::Apply() {

  // default methods to be trained + tested
  cout << "    Starting TMVA application:" << endl;

  // create reader
  Reader* reader = new Reader(_sReadOpt.data());
  for (const string variable : _vecVarsForTMVA) {
    reader -> AddVariable(variable.data(), &_mapInTupleVars[variable]);
  }
  cout << "      Created reader." << endl;

  // book methods to use
  for (const auto methodAndOpts : _vecMethodsAndOptsTMVA) {
    const string name    = get<0>(methodAndOpts);
    const string weights = _sLoader + "/weights/" + _sFactory + "_" + get<0>(methodAndOpts) + ".weights.xml";
    reader -> BookMVA(name, weights);
  }
  cout << "      Booked methods." << endl;

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
      // TODO extend to multiple targets
      const float target = (reader -> EvaluateRegression(get<0>(methodAndOpt)))[0];
      _mapTmvaHists[get<0>(methodAndOpt)] -> Fill(target);

      /* TODO fill output histograms */

      // form leaf and it to list
      const string sLeaf = get<0>(methodAndOpt) + "_" + _vecTargsForTMVA[0];

      // initialize value of leaf
      _mapOutTupleVars[sLeaf] = target;;
    }  // end method loop

    // fill output tuple with regression targets & input leaves
    for (const auto inBranch : _mapInTupleVars) {
      _mapOutTupleVars[inBranch.first] = inBranch.second;
    }
    FillTuples();

  }  // end event loop
  cout << "    Application finished!" << endl;

  // delete reader and exit
  delete reader;
  return;

}  // end 'Apply()'



void BHCalCalibration::End() {

  // announce finishing
  cout << "    Finishing:" << endl;

  FillHistos();
  ComputeReso();
  SaveOutput();

  cout << "  Finished BHCal calibration!\n" << endl;
  return;

}  // end 'End()'



// setters --------------------------------------------------------------------

void BHCalCalibration::SetInput(const string sInFile, const string sInTuple, const float treeWeight) {

  _sInFile  = sInFile;
  _sInTuple = sInTuple;
  _weight   = treeWeight;
  return;

}  // end 'SetInput(string, string, float)'



void BHCalCalibration::SetTupleLeaves(const vector<string> vecInput) {

  _vecInTupleLeaves = vecInput;
  for (const string inputLeaf : _vecInTupleLeaves) {
    _mapInTupleVars[inputLeaf] = -999.;
  }
  return;

}  // end 'SetTupleArgs(vector<string>)'



void BHCalCalibration::SetTmvaOpts(const string sFactOpt, const string sTrainOpt, const string sReadOpt, const bool addSpecs) {

  _sFactOpt  = sFactOpt;
  _sTrainOpt = sTrainOpt;
  _sReadOpt  = sReadOpt;
  _addSpecs  = addSpecs;
  return;

}  // end 'SetTmvaOpts(string, string, string, bool)'



void BHCalCalibration::SetTmvaArgs(const vector<string> vecVars, const vector<string> vecTargs, const vector<string> vecSpecs, const TCut select) {

  _vecVarsForTMVA  = vecVars;
  _vecTargsForTMVA = vecTargs;
  _vecSpecsForTMVA = vecSpecs;
  _cSelect         = select;
  return;

}  // end 'SetTmvaArgs(vector<string>, vector<string>, vector<string>, TCut)'



void BHCalCalibration::SetTmvaMethods(const vector<pair<string, string>> vecMethodAndOpts) {

  for (const auto methodAndOpt : vecMethodAndOpts) {
    _vecMethodsAndOptsTMVA.push_back(make_tuple(methodAndOpt.first, methodAndOpt.second, _mapMethodToIndex[methodAndOpt.first]));
  }

  /* TODO check for duplicate methods */
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

  cout << "      Opened files." << endl;
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

  // set input branches and collect leaves into a list
  for (const string leaf : _vecInTupleLeaves) {
    _ntInput -> SetBranchAddress(leaf.data(), &_mapInTupleVars[leaf]);
  }

  // add targets to output leaf list
  string sLeaf     = "";
  string sLeafList = "";
  for (const auto methodAndOpt : _vecMethodsAndOptsTMVA) {
    for (const string target : _vecTargsForTMVA) {

      // form leaf and add it to list
      sLeaf      = target + "_" + get<0>(methodAndOpt);
      sLeafList += sLeaf;
      sLeafList += ":";

      // initialize value of leaf
      _mapOutTupleVars[sLeaf] = -999.;
    }
  }

  // determine total no. of input leaves
  const size_t nLeaves = _vecInTupleLeaves.size();

  // add input leaves to output list
  size_t iLeaf = 0;
  for (const string leaf : _vecInTupleLeaves) {

    // form leaf and add it to list
    sLeafList += leaf;
    if ((iLeaf + 1) != nLeaves) {
      sLeafList += ":";
    }

    _mapOutTupleVars[leaf] = -999.;
    ++iLeaf;
  }

  _ntOutput = new TNtuple(_sOutTuple.data(), "regression targets vs. input", sLeafList.data());
  cout << "      Initialized tuples." << endl;
  return;

}  // end 'InitTuples()'



void BHCalCalibration::FillTuples() {

  // load regression targets into a vector
  _vecOutTupleValues.clear();
  for (const auto outLeaf : _mapOutTupleVars) {
    _vecOutTupleValues.push_back(outLeaf.second);
  }

  // fill output tuple
  _ntOutput -> Fill(_vecOutTupleValues.data());
  return;

}  // end 'FillTuples()'



void BHCalCalibration::InitHistos() {

  // binning for tmva book keeping histograms
  const tuple<size_t, float, float> tmvaBins = make_tuple(100, -100., 600.);

  // create tmva book keeping histograms
  // TODO extend to multiple targets
  for (const auto methodAndOpt : _vecMethodsAndOptsTMVA) {
    const string sTmvaName = "h_" + get<0>(methodAndOpt);
    _mapTmvaHists[get<0>(methodAndOpt)] = new TH1F(sTmvaName.data(), get<0>(methodAndOpt).data(), get<0>(tmvaBins), get<1>(tmvaBins), get<2>(tmvaBins));
  }

  cout << "      Initialized histograms." << endl;
  return;

}  // end 'InitHistos()'



void BHCalCalibration::FillHistos() {

  /* TODO add histogram filling */
  cout << "      Filled histograms." << endl;
  return;

}  // end 'FillHistos()'



void BHCalCalibration::ComputeReso() {

  /* TODO add resolution computation */
  cout << "      Calculated resolutions." << endl;
  return;

}  // end 'ComputeReso()'



void BHCalCalibration::SaveOutput() {

  // create directories
  TDirectory *dTMVA = (TDirectory*) _fOutput -> mkdir("tmva");

  // save tmva histograms
  dTMVA -> cd();
  for (auto tmvaHist : _mapTmvaHists) {
    tmvaHist.second -> Write();
  }

  // save output tuple
  _fOutput  -> cd();
  _ntOutput -> Write();

  cout << "      Saved histograms." << endl;
  return;

}  // end 'SaveOutput()'



void BHCalCalibration::CloseFiles() {

  _fOutput -> cd();
  _fOutput -> Close();
  _fInput  -> cd();
  _fInput  -> Close();

  cout << "     Closed files." << endl;
  return;

}  // end 'CloseFiles()'

// end ------------------------------------------------------------------------
