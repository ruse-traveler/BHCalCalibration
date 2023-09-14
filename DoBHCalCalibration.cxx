// ----------------------------------------------------------------------------
// 'DoBHCalCalibration.cc'
// Derek Anderson
// 09.14.2023
//
// Driver macro to run ePIC BHCal calibration workflow.
// ----------------------------------------------------------------------------

// c utilities
#include <string>
#include <vector>
#include <utility>
// root utilites
#include <TROOT.h>
#include <TSystem.h>
// calibration class
#include "src/BHCalCalibration.hh"

// load libraries
R__LOAD_LIBRARY(./src/BHCalCalibration_cc.so)

using namespace std;



void DoBHCalCalibration() {

  // load library
  //gROOT -> ProcessLine(".L ./src/BHCalCalibration.cc++");
  //gSystem -> Load("./src/BHCalCalibration_cc.so");

  /* TODO options go here */

  // run calibration workflow
  BHCalCalibration* calibrator = new BHCalCalibration();
  calibrator -> Init();
  calibrator -> Train();
  calibrator -> Apply();
  calibrator -> End();
  return;

}

// end ------------------------------------------------------------------------
