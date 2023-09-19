// ----------------------------------------------------------------------------
// 'BHCalCalibration.hh'
// Derek Anderson
// 09.14.2023
//
// A simple interface to calibrate the simulated ePIC BHCal response with
// TMVA and calculate its energy resolution.
// ----------------------------------------------------------------------------

#ifndef BHCALCALIBRATION_HH
#define BHCALCALIBRATION_HH

// c utilities
#include <map>
#include <vector>
#include <string>
#include <utility>
// root utilities
#include <TCut.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TBranch.h>
// tmva classes
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"

using namespace std;



// BHCalCalibration definition ------------------------------------------------

class BHCalCalibration {

  public:

    // ctor/dtor
    BHCalCalibration(const string sFactory = "TMVARegression", const string sLoader = "TMVADir", const string sOutput = "out.root")  {};
    ~BHCalCalibration();

    // analysis methods
    void Init();
    void Train();
    void Apply();
    void End();

    // setters
    void SetInput(const string sInFile, const string sInTuple, const float treeWeight = 1.);
    void SetTupleArgs(const vector<string> vecInput);
    void SetTmvaOpts(const string sFactOpt, const string sTrainOpt, const string sReadOpt, const bool addSpecs = false);
    void SetTmvaArgs(const vector<string> vecVars, const vector<string> vecTargs, const vector<string> vecSpecs = {}, const TCut select = "");
    void SetTmvaMethods(const vector<pair<string, string>> vecMethodAndOpts);

  private:

    // private methods
    void OpenFiles();
    void InitTuples();
    void InitHistos();
    void FillHistos();
    void ComputeReso();
    void SaveOutput();
    void CloseFiles();

    // i/o members
    string   _sOutFile = "";
    string   _sInFile  = "";
    string   _sInTuple = "";
    TFile*   _fInput   = NULL;
    TFile*   _fOutput  = NULL;
    TNtuple* _ntInput  = NULL;
    TNtuple* _ntOutput = NULL;

    // tuple members
    vector<string>     _vecInTupleLeaves;
    vector<string>     _vecOutTupleLeaves;
    map<string, float> _mapInTupleVars;
    map<string, float> _mapOutTupleVars;

    // general tmva parameters
    bool   _addSpecs  = false;
    float  _weight    = 1.;
    string _sFactory  = "TMVARegression";
    string _sLoader   = "TMVADir";
    string _sFactOpt  = "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression";
    string _sTrainOpt = "nTrain_Regression=1000:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V";
    string _sReadOpt  = "!Color:!Silent";

    // tmva training arguments
    TCut           _cSelect = "";
    vector<string> _vecVarsForTMVA;
    vector<string> _vecTargsForTMVA;
    vector<string> _vecSpecsForTMVA;

    // tmva methods, options, and types
    vector<tuple<string, string, TMVA::Types>> _vecMethodsAndOptsTMVA;

    // method-index map
    map<string, TMVA::Types> _mapMethodToIndex;

};  // end BHCalCalibration definition



// ctor/dtor ------------------------------------------------------------------

BHCalCalibration::BHCalCalibration(const string sFactory, const string sLoader, const string sOutput) {

  _sFactory = sFactory;
  _sLoader  = sLoader;
  _sOutFile = sOutput;

  // make sure vectors are clear
  _vecInTupleLeaves.clear();
  _vecOutTupleLeaves.clear();
  _vecVarsForTMVA.clear();
  _vecTargsForTMVA.clear();
  _vecSpecsForTMVA.clear();
  _vecMethodsAndOptsTMVA.clear();

  // make sure maps are clear
  _mapInTupleVars.clear();
  _mapOutTupleVars.clear();

  // create method-index map
  _mapMethodToIndex = {
    {"PDERS",    TMVA::Types::kPDERS},
    {"PDEFoam",  TMVA::Types::kPDEFoam},
    {"KNN",      TMVA::Types::kKNN},
    {"LD",       TMVA::Types::kLD},
    {"FDA_MC",   TMVA::Types::kFDA},
    {"FDA_GA",   TMVA::Types::kFDA},
    {"FDA_GAMT", TMVA::Types::kFDA},
    {"MLP",      TMVA::Types::kMLP},
    {"DNN_CPU",  TMVA::Types::kDL},
    {"DNN_GPU",  TMVA::Types::kDL},
    {"SVM",      TMVA::Types::kSVM},
    {"BDT",      TMVA::Types::kBDT},
    {"BDTG",     TMVA::Types::kBDT}
  };
  cout << "CTOR" << endl;

}  // end ctor



BHCalCalibration::~BHCalCalibration() {

  delete _fInput;
  delete _fOutput;
  delete _ntInput;
  delete _ntOutput;
  cout << "DTOR" << endl;

}  // end dtor

#endif

// end ------------------------------------------------------------------------
