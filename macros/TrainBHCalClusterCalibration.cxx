/// ===========================================================================
/*! \file   TrainBHCalClusterCalibration.cxx
 *  \author Derek Anderson
 *  \date   09.17.2024
 *
 *  A ROOT macro to train a TMVA model to calibrate
 *  the energy of clusters in the BHCal and BIC.
 *  Ingests the TNtuple produced by
 *  'FillBHCalClusterCalibrationTuple.cxx'
 */
/// ===========================================================================

#define TrainBHCalClusterCalibration_cxx

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
// analysis utilities
#include "../../utility/TMVAHelper.hxx"



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
  "testA.root",
  "tmva_test",
  "TMVARegression",
  true
};



// ============================================================================
//! Train a TMVA model for BHCal cluster calibration
// ============================================================================
void TrainBHCalClusterCalibration(const Options& opt = DefaultOptions) {

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

  // other tmva options
  const bool  addSpectators(false);
  const bool  doECalCut(false);
  const float treeWeight(1.0);
  const TCut  trainCut("(eSumBHCal>=0)&&(eSumBEMC>=0)&&(abs(hLeadBHCal)<1.1)&&(abs(hLeadBEMC)<1.1)");

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
  if (!ntInput) {
    std::cerr << "PANIC: couldn't grab input tuple!\n"
              << "       name  = " << opt.in_tuple << "\n"
              << "       tuple = " << ntInput
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
  TMVAHelper::Trainer train_helper(vecUseAndVar, vecMethods);
  train_helper.SetFactoryOptions(vecFactoryOpts);
  train_helper.SetTrainOptions(vecTrainOpts);
  std::cout << "    Created TMVA helper." << std::endl;


  // --------------------------------------------------------------------------
  // Train tmva models
  // --------------------------------------------------------------------------

  // instantiate tmva library
  TMVA::Tools::Instance();
  std::cout << "    Begin training calibration models:" << std::endl;

  // create tmva factory & load data
  TMVA::Factory*    factory = new TMVA::Factory(opt.name_tmva.data(), output, train_helper.CompressFactoryOptions().data());
  TMVA::DataLoader* loader  = new TMVA::DataLoader(opt.out_tmva.data());
  std::cout << "      Created factory and data loader..." << std::endl;

  // now load variables
  train_helper.LoadVariables(loader, addSpectators);
  std::cout << "      Loaded variables..." << std::endl;

  // add tree & prepare for training
  loader -> AddRegressionTree(ntInput, treeWeight);
  loader -> PrepareTrainingAndTestTree(trainCut, train_helper.CompressTrainingOptions().data());
  std::cout << "      Added tree, prepared training..." << std::endl;

  // book methods
  train_helper.BookMethodsToTrain(factory, loader);
  std::cout << "      Booked methods for training..." << std::endl;

  // train, test, & evaluate
  factory -> TrainAllMethods();
  factory -> TestAllMethods();
  factory -> EvaluateAllMethods();
  std::cout << "      Trained models.\n"
            << "    Finished training calibration models!"
            << endl;

  // --------------------------------------------------------------------------
  // Close I/O and exit
  // --------------------------------------------------------------------------

  // close files
  output -> cd();
  output -> Close();
  input  -> cd();
  input  -> Close();

  // delete tmva objects
  delete factory;
  delete loader;

  // announce end & exit
  std::cout << "  Finished BHCal training script!\n" << std::endl;
  return;

}

// end ========================================================================
