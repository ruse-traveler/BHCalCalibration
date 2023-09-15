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
#include <vector>
#include <string>
#include <utility>
// root utilities
#include <TCut.h>
#include <TFile.h>
#include <TNtuple.h>

using namespace std;


// BHCalCalibration definition ------------------------------------------------

class BHCalCalibration {

  public:

    // ctor/dtor
    BHCalCalibration()  {};
    ~BHCalCalibration() {};

    // public methods
    void Init();
    void Train();
    void Apply();
    void End();

    // setters
    // TODO finish adding TMVA setters
    void SetTupleArgs(vector<string> vecInput);
    void SetTmvaOpts(string sFactOpt, string sReadOpt, string 
    void SetTmvaArgs(vector<string> vecVars, vector<string> vecSpecs, vector<string> vecTars);
    void SetTmvaMethdods(pair<vector<string>, vector<string>> vecMethodAndOpts);

  private:

    // private methods
    void ParseInput();
    void InitTuples();
    void InitHistos();
    void FillHistos();
    void ComputeReso();
    void SaveOutput();

    // i/o members
    TFile*   _fInput   = NULL;
    TFile*   _fOutput  = NULL;
    TNtuple* _ntInput  = NULL;
    TNtuple* _ntOutput = NULL;

    // tuple members
    vector<float>  _vecInTupleVars;
    vector<float>  _vecOutTupleVars;
    vector<string> _vecInTupleLeaves;
    vector<string> _vecOutTupleLeaves;

    // tmva parameters
    // TODO add TMVA Type enum to input methods
    TCut           _cSelect   = "";
    bool           _addSpecs  = false;
    float          _weight    = 1.;
    string         _sFactory  = "TMVARegression";
    string         _sFactOpt  = "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression";
    string         _sTrainOpt = "nTrain_Regression=1000:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V";
    string         _sReadOpt  = "!Color:!Silent";
    string         _sLoader   = "";
    vector<string> _vecVarsForTMVA;
    vector<string> _vecSpecForTMVA;
    vector<string> _vecTarsForTMVA;
    vector<string> _vecMethodsTMVA;
    vector<string> _vecMethodOptsTMVA;

};  // end BHCalCalibration definition

#endif

// end ------------------------------------------------------------------------
