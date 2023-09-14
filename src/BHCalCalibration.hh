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
#include <utility>
// root utilities
#include <TFile.h>
#include <TNtuple.h>


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

};  // end BHCalCalibration definition

#endif

// end ------------------------------------------------------------------------
