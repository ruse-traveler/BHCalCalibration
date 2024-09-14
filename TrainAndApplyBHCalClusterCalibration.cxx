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
  std::string out_tmva;     // output tmva directory
  std::string in_file;      // input file
  std::string in_tuple;     // input ntuple
  bool        do_progress;  // print progress through entry loop
}  DefaultOptions = {
  "test_v0.root",
  "tmva_test",
  "./input/forRefinedTmvaMacro_calibrationTupleWithNEvt1K.d12m9y2024.root",
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

  // methods to use
  //   - TODO might be good to add field for method-specific options
  const std::vector<std::string> vecMethods = {
    "LD",
    "kNN",
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
    "nTrain_Regression=10",
    "nTest_Regression=0",
    "SplitMode=Random:NormMode=NumEvents",
    "!V"
  };
  const std::vector<std::string> vecReadOpts = {
    "!Color",
    "!Silent"
  };

  // other tmva options
  //   - TODO clean up
  const bool  addSpectators(false);
  const float treeWeight(1.0);
  const TCut  trainCut("(eSumBHCal>0)&&(eSumBEMC>=0)&&(abs(hLeadBHCal)<1.1)&&(abs(hLeadBEMC)<1.1)");

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

  // --------------------------------------------------------------------------
  // Set up helpers
  // --------------------------------------------------------------------------

  // create tmva helper
  TMVAHelper tmva_helper(vecUseAndVar, vecMethods);
  tmva_helper.SetFactoryOptions(vecFactoryOpts);
  tmva_helper.SetTrainOptions(vecTrainOpts);
  tmva_helper.SetReadOptions(vecReadOpts);

  // lambda to collect input leaves into a single vector
  auto getLeaves = [](const std::vector<std::pair<TMVAHelper::Use, std::string>>& inputs) {
    std::vector<std::string> outputs;
    for (const auto& input : inputs) {
      outputs.push_back(input.second);
    }
    return outputs;
  };

  // create input/output helpers
  NTupleHelper in_helper( getLeaves(vecUseAndVar) );
  NTupleHelper out_helper( tmva_helper.GetOutputs() );

  // set input/output tuple branches
  TNtuple* ntOutput = new TNtuple("ntTmvaOutput", "Output of TMVA regression", out_helper.CompressVariables().data());
  in_helper.SetBranches(ntToApply);
  std::cout << "    Set input/output tuple branches." << std::endl;

  // --------------------------------------------------------------------------
  // Train tmva models
  // --------------------------------------------------------------------------

  // instantiate tmva library
  TMVA::Factory*    factory;
  TMVA::DataLoader* loader;
  TMVA::Tools::Instance();
  std::cout << "    Begin training calibration models:" << std::endl;

  // create tmva factory & load data
  factory = new TMVA::Factory("TMVARegression", output, tmva_helper.CompressFactoryOptions().data());
  loader  = new TMVA::DataLoader(opt.out_tmva.data());
  std::cout << "      Created factory and data loader..." << std::endl;

  // now load variables
  tmva_helper.LoadVariables(loader, addSpectators);
  std::cout << "      Loaded variables..." << std::endl;

  // add tree & prepare for training
  loader -> AddRegressionTree(ntToTrain, treeWeight);
  loader -> PrepareTrainingAndTestTree(trainCut, tmva_helper.CompressTrainOptions().data());
  std::cout << "      Added tree, prepared training..." << std::endl;

  // book methods
  tmva_helper.BookMethodsToTrain(factory, loader);
  std::cout << "      Booked methods for training..." << std::endl;

  // train, test, & evaluate
  factory -> TrainAllMethods();
  factory -> TestAllMethods();
  factory -> EvaluateAllMethods();
  std::cout << "      Trained models.\n"
            << "    Finished training calibration models!"
            << endl;

  // --------------------------------------------------------------------------
  // Apply tmva models
  // --------------------------------------------------------------------------

  // instantiate reader
  TMVA::Reader* reader = new Reader(tmva_helper.CompressReadOptions().data());
  std::cout << "    Begin applying calibration models:" << std::endl;

  // add input variables to reader, book methods
  tmva_helper.AddVariables(reader, in_helper);
  tmva_helper.BookMethodsToRead(reader, opt.out_tmva);
  std::cout << "      Added variables and methods to read." << std::endl;

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

    // make sure output variables are empty
    out_helper.ResetValues();

    /* TODO evaluate targets, fill output tuple */

  }  // end entry loop
  std::cout << "    Application loop finished." << std::endl;

  // close files
  output -> cd();
  output -> Close();
  input  -> cd();
  input  -> Close();

  // delete tmva objects
  delete factory;
  delete loader;
  delete reader;

  // announce end & exit
  std::cout << "  Finished BHCal calibration script!\n" << std::endl;
  return;

}

// end ========================================================================
