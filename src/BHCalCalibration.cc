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
#include <iostream>
// class header
#include "BHCalCalibration.hh"

using namespace std;



// public methods -------------------------------------------------------------

void BHCalCalibration::Init() {

  ParseInput();
  InitTuples();
  InitHistos();

  cout << "INIT" << endl;
  return;

}  // end 'Init()'



void BHCalCalibration::Train() {

  cout << "TRAIN" << endl;
  return;

}  // end 'Train()'



void BHCalCalibration::Apply() {

  cout << "APPLY" << endl;
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
