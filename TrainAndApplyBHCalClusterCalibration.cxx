/// ===========================================================================
/*! \file   TrainAndApplyBHCalClusterCalibration.cxx
 *  \author Derek Anderson
 *  \date   09.11.2024
 *
 *  A ROOT macro to train and apply a TMVA model to
 *  calibrate the energy of clusters in the BHCal
 *  and BIC. Ingests the TNtuple produced by
 *  'FillBHCalClusterCalibrationTuple.cxx'
 */
/// ===========================================================================

#define TrainAndApplyBHCalClusterCalibration_cxx

// c++ utilities
#include <string>
#include <vector>
#include <cassert>
#include <utility>
#include <iostream>
// root libraries
#include <TCut.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TSystem.h>
// tmva components
#include <TMVA/Tools.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/TMVARegGui.h>
// analysis utilities
#include "../../utility/TMVAHelper.hxx"
#include "../../utility/NTupleHelper.hxx"



// ============================================================================
//! Struct to consolidate user options
// ============================================================================
struct Options {
  std::string out_file;     // output file
  std::string in_file;      // input file
  std::string in_tuple;     // input ntuple
  bool        do_progress;  // print progress through entry loop
}  DefaultOptions = {
  "test_v0.root",
  "./input/testRefinedTmvaMacro_calibrationTuple.d10m9y2024.root",
  "ntForCalib",
  true
};



// ============================================================================
//! Train and apply a TMVA model for BHCal cluster calibration
// ============================================================================
void TrainAndApplyBHCalClusterCalibration(const Options& opt = DefaultOptions) {

  // --------------------------------------------------------------------------
  // calculation parameters
  // --------------------------------------------------------------------------

  // input variables & usage
  std::vector<std::pair<TMVAHelper::Use, std::string>> vecUseAndVar = {
    {TMVAHelper::Use::Target, "ePar"},
    {TMVAHelper::Use::Watch, "fracParVsLeadBHCal"},
    {TMVAHelper::Use::Watch, "fracParVsLeadBEMC"},
    {TMVAHelper::Use::Watch, "fracParVsSumBHCal"},
    {TMVAHelper::Use::Watch, "fracParVsSumBEMC"},
    {TMVAHelper::Use::Watch, "fracLeadBHCalVsBEMC"},
    {TMVAHelper::Use::Watch, "fracSumBHCalVsBEMC"},
    {TMVAHelper::Use::Train, "eLeadBHCal"},
    {TMVAHelper::Use::Train, "eLeadBEMC"},
    {TMVAHelper::Use::Watch, "eSumBHCal"},
    {TMVAHelper::Use::Watch, "eSumBEMC"},
    {TMVAHelper::Use::Watch, "diffLeadBHCal"},
    {TMVAHelper::Use::Watch, "diffLeadBEMC"},
    {TMVAHelper::Use::Watch, "diffSumBHCal"},
    {TMVAHelper::Use::Watch, "diffSumBEMC"},
    {TMVAHelper::Use::Train, "nHitsLeadBHCal"},
    {TMVAHelper::Use::Train, "nHitsLeadBEMC"},
    {TMVAHelper::Use::Watch, "nClustBHCal"},
    {TMVAHelper::Use::Watch, "nClustBEMC"},
    {TMVAHelper::Use::Train, "hLeadBHCal"},
    {TMVAHelper::Use::Train, "hLeadBEMC"},
    {TMVAHelper::Use::Train, "fLeadBHCal"},
    {TMVAHelper::Use::Train, "fLeadBEMC"},
    {TMVAHelper::Use::Watch, "eLeadImage"},
    {TMVAHelper::Use::Train, "eSumImage"},
    {TMVAHelper::Use::Watch, "eLeadScFi"},
    {TMVAHelper::Use::Train, "eSumScFi"},
    {TMVAHelper::Use::Watch, "nClustImage"},
    {TMVAHelper::Use::Watch, "nClustScFi"},
    {TMVAHelper::Use::Watch, "hLeadImage"},
    {TMVAHelper::Use::Watch, "hLeadScFi"},
    {TMVAHelper::Use::Watch, "fLeadImage"},
    {TMVAHelper::Use::Watch, "fLeadScFi"},
    {TMVAHelper::Use::Train, "eSumScFiLayer1"},
    {TMVAHelper::Use::Train, "eSumScFiLayer2"},
    {TMVAHelper::Use::Train, "eSumScFiLayer3"},
    {TMVAHelper::Use::Train, "eSumScFiLayer4"},
    {TMVAHelper::Use::Train, "eSumScFiLayer5"},
    {TMVAHelper::Use::Train, "eSumScFiLayer6"},
    {TMVAHelper::Use::Train, "eSumScFiLayer7"},
    {TMVAHelper::Use::Train, "eSumScFiLayer8"},
    {TMVAHelper::Use::Train, "eSumScFiLayer9"},
    {TMVAHelper::Use::Train, "eSumScFiLayer10"},
    {TMVAHelper::Use::Train, "eSumScFiLayer11"},
    {TMVAHelper::Use::Train, "eSumScFiLayer12"},
    {TMVAHelper::Use::Train, "eSumImageLayer1"},
    {TMVAHelper::Use::Watch, "eSumImageLayer2"},
    {TMVAHelper::Use::Train, "eSumImageLayer3"},
    {TMVAHelper::Use::Train, "eSumImageLayer4"},
    {TMVAHelper::Use::Watch, "eSumImageLayer5"},
    {TMVAHelper::Use::Train, "eSumImageLayer6"}
  };

  /* TODO add remaining parameters */

  // lower verbosity & announce start
  gErrorIgnoreLevel = kError;
  std::cout << "\n  Beginning calibration training and evaluation macro..." << std::endl;

  // --------------------------------------------------------------------------
  // Open input/outputs
  // --------------------------------------------------------------------------

  // open files
  TFile* input  = new TFile(opt.in_file.data(),  "read");
  TFile* output = new TFile(opt.out_file.data(), "recreate");
  if (!input || !output) {
    std::cerr << "PANIC: couldn't open a file!\n"
              << "       input  = " << input << "\n"
              << "       output = " << output
              << std::endl;
    assert(input && output);
  }

  // print input/output files
  std::cout << "    Opened input/output files:\n"
            << "      input file  = " << opt.in_file << "\n"
            << "      output file = " << opt.out_file
            << std::endl;

  // grab input tuple
  TNtuple* ntToTrain = (TNtuple*) input -> Get(opt.in_tuple.data());
  TNtuple* ntToApply = (TNtuple*) input -> Get(opt.in_tuple.data());
  if (!ntToTrain || !ntToApply) {
    std::cerr << "PANIC: couldn't grab input tuple!\n"
              << "       name  = " << opt.in_tuple << "\n"
              << "       train = " << ntToTrain    << "\n"
              << "       apply = " << ntToApply
              << std::endl;
    assert(ntToTrain && ntToApply);
  }

  // print input tuples
  std::cout << "    Grabbed input tuples:\n"
            << "      tuple = " << opt.in_tuple
            << std::endl;

  // create tmva helper
  //   - TODO fill in remaining setters/etc
  //   - TODO generate list ouf output leaves
  TMVAHelper tmva_helper;
  tmva_helper.SetInputVariables(vecUseAndVar);

  // create input/output helpers
  //   - TODO fill in
  NTupleHelper in_helper;
  NTupleHelper out_helper;

  // set input tuple branches
  in_helper.SetBranches(ntToApply);
  std::cout << "    Set tuple branches." << std::endl;

  // --------------------------------------------------------------------------
  // Train tmva models
  // --------------------------------------------------------------------------

  /* TODO fill in */

  // --------------------------------------------------------------------------
  // Apply tmva models
  // --------------------------------------------------------------------------

  /* TODO set up reader, etc. */

  // get number of events for application
  const uint64_t nEntries = ntToApply -> GetEntries();
  cout << "    Processing: " << nEntries << " events" << endl;

  // loop over input tuple
  uint64_t nBytes = 0;
  for (uint64_t iEntry = 0; iEntry < nEntries; iEntry++) {

    // announce progress
    if (opt.do_progress) {
      std::cout << "      Processing entry " << iEntry + 1 << "/" << nEntries << "...";
      if (iEntry + 1 < nEntries) {
        std::cout << "\r" << std::flush;
      } else {
        std::cout << std::endl;
      }
    }

    // grab entry
    const uint64_t bytes = ntToApply -> GetEntry(iEntry);
    if (bytes < 0.) {
      std::cerr << "WARNING error in entry #" << iEntry << "! Aborting loop!" << std::endl;
      break;
    } else {
      nBytes += bytes;
    }

    /* TODO evaluate targets, fill output tuple */

  }  // end entry loop
  std::cout << "    Application loop finished." << std::endl;

  // close files
  output -> cd();
  output -> Close();
  input  -> cd();
  input  -> Close();

  // announce end & exit
  std::cout << "  Finished BHCal calibration script!\n" << std::endl;
  return;

}

// end ========================================================================
