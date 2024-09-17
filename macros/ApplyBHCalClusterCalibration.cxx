/// ===========================================================================
/*! \file   ApplyBHCalClusterCalibration.cxx
 *  \author Derek Anderson
 *  \date   09.17.2024
 *
 *  A ROOT macro to apply a pre-trained TMVA model to
 *  calibrate the energy of clusters in the BHCal and
 *  BIC. Ingests the TNtuple produced by
 *  'FillBHCalClusterCalibrationTuple.cxx'
 */
/// ===========================================================================

#define ApplyBHCalClusterCalibration_cxx

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
#include <TMVA/Reader.h>
// analysis utilities
#include "../../utility/TMVAHelper.hxx"
#include "../../utility/NTupleHelper.hxx"



// ============================================================================
//! Struct to consolidate user options
// ============================================================================
struct Options {
  std::string in_file;      // input file
  std::string in_tuple;     // input ntuple
  std::string out_file;     // output file
  std::string out_tmva;     // output tmva directory
  std::string name_tmva;    // name of TMVA process
  bool        do_progress;  // print progress through entry loop
}  DefaultOptions = {
  "./input/forNewTrainingMacro_noNonzeroEvts_andDefinitePrimary.evt5Ke210pim_central.d14m9y2024.root",
  "ntForCalib",
  "testB.root",
  "tmva_test",
  "TMVARegression",
  true
};



// ============================================================================
//! Apply a TMVA model for BHCal cluster calibration
// ============================================================================
void ApplyBHCalClusterCalibration(const Options& opt = DefaultOptions) {

  // --------------------------------------------------------------------------
  // calculation parameters
  // --------------------------------------------------------------------------

  // input variables & usage
  const std::vector<std::pair<TMVAHelper::Use, std::string>> vecUseAndVar = {
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
    {TMVAHelper::Use::Watch, "nHitsLeadBHCal"},
    {TMVAHelper::Use::Watch, "nHitsLeadBEMC"},
    {TMVAHelper::Use::Watch, "nClustBHCal"},
    {TMVAHelper::Use::Watch, "nClustBEMC"},
    {TMVAHelper::Use::Watch, "hLeadBHCal"},
    {TMVAHelper::Use::Watch, "hLeadBEMC"},
    {TMVAHelper::Use::Watch, "fLeadBHCal"},
    {TMVAHelper::Use::Watch, "fLeadBEMC"},
    {TMVAHelper::Use::Watch, "eLeadImage"},
    {TMVAHelper::Use::Watch, "eSumImage"},
    {TMVAHelper::Use::Watch, "eLeadScFi"},
    {TMVAHelper::Use::Watch, "eSumScFi"},
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

  // methods to use
  //   - TODO might be good to add field for method-specific options
  //     (e.g. FDA_GA needs a function to work)
  const std::vector<std::string> vecMethods = {
    "LD",
    "KNN",
    "MLP",
    "BDTG",
    "FDA_GA",
    "PDEFoam"
  };

  // general tmva options
  const std::vector<std::string> vecFactoryOpts = {
    "!V",
    "!Silent",
    "Color",
    "DrawProgressBar",
    "AnalysisType=Regression"
  };
  const std::vector<std::string> vecTrainOpts = {
    "nTrain_Regression=100",
    "nTest_Regression=0",
    "SplitMode=Random:NormMode=NumEvents",
    "!V"
  };
  const std::vector<std::string> vecReadOpts = {
    "!Color",
    "!Silent"
  };

  // other tmva options
  const bool  addSpectators(false);
  const bool  doECalCut(false);
  const float treeWeight(1.0);
  const TCut  trainCut("(eSumBHCal>=0)&&(eSumBEMC>=0)&&(abs(hLeadBHCal)<1.1)&&(abs(hLeadBEMC)<1.1)");

  // ecal cut parameters
  const std::pair<double, double> eneECalRange = {0.5, 100.};

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
    assert(output && input);
  }

  // print input/output files
  std::cout << "    Opened input/output files:\n"
            << "      input file  = " << opt.in_file << "\n"
            << "      output file = " << opt.out_file
            << std::endl;

  // grab input tuple
  TNtuple* ntInput = (TNtuple*) input -> Get(opt.in_tuple.data());
  if (!input) {
    std::cerr << "PANIC: couldn't grab input tuple!\n"
              << "       name  = " << opt.in_tuple << "\n"
              << "       input = " << input
              << std::endl;
    assert(ntInput);
  }

  // print input tuples
  std::cout << "    Grabbed input tuples:\n"
            << "      tuple = " << opt.in_tuple
            << std::endl;

  // --------------------------------------------------------------------------
  // Set up helpers
  // --------------------------------------------------------------------------

  // create tmva helper
  TMVAHelper::Reader read_helper(vecUseAndVar, vecMethods);
  read_helper.SetOptions(vecReadOpts);

  // collect input leaves into a single vector
  std::vector<std::string> inputs;
  for (const auto& useAndVar : vecUseAndVar) {
    inputs.push_back(useAndVar.second);
  }

  // create input/output helpers
  NTupleHelper in_helper( inputs );
  NTupleHelper out_helper( read_helper.GetOutputs() );

  // set input/output tuple branches
  TNtuple* ntOutput = new TNtuple("ntTmvaOutput", "Output of TMVA regression", out_helper.CompressVariables().data());
  in_helper.SetBranches(ntInput);
  std::cout << "    Set input/output tuple branches." << std::endl;

  // --------------------------------------------------------------------------
  // Apply tmva models
  // --------------------------------------------------------------------------

  // instantiate reader
  TMVA::Reader* reader = new TMVA::Reader(read_helper.CompressOptions().data());
  std::cout << "    Begin applying calibration models:" << std::endl;

  // add input variables to reader, book methods
  read_helper.ReadVariables(reader, in_helper);
  read_helper.BookMethodsToRead(reader, opt.out_tmva, opt.name_tmva);
  std::cout << "      Added variables and methods to read." << std::endl;

  // get number of events for application
  const uint64_t nEntries = ntInput -> GetEntries();
  cout << "    Processing: " << nEntries << " events" << endl;

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
    const uint64_t bytes = ntInput -> GetEntry(iEntry);
    if (bytes < 0.) {
      std::cerr << "WARNING error in entry #" << iEntry << "! Aborting loop!" << std::endl;
      break;
    } else {
      nBytes += bytes;
    }

    // make sure output variables are empty
    out_helper.ResetValues();
    read_helper.ResetValues();

    // evaluate targets
    read_helper.EvaluateMethods(reader, in_helper);

    // apply ecal cut if need be
    const double eLeadBEMC   = in_helper.GetVariable("eLeadBEMC");
    const bool   isInECalCut = ((eLeadBEMC > eneECalRange.first) && (eLeadBEMC < eneECalRange.second));
    if (doECalCut && !isInECalCut) continue;

    // set values in output tuple & fill
    for (const std::string& output : read_helper.GetOutputs()) {
      out_helper.SetVariable( output, read_helper.GetVariable(output) );
    }
    ntOutput -> Fill( out_helper.GetValues().data() );

  }  // end entry loop
  std::cout << "    Application loop finished." << std::endl;

  // --------------------------------------------------------------------------
  // Save output and exit
  // --------------------------------------------------------------------------

  // save & close files
  output   -> cd();
  ntOutput -> Write(); 
  output   -> Close();
  input    -> cd();
  input    -> Close();

  // delete tmva object
  delete reader;

  // announce end & exit
  std::cout << "  Finished BHCal calibration script!\n" << std::endl;
  return;

}

// end ========================================================================
