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
#include <TTreeFormula.h>
// tmva components
#include <TMVA/Tools.h>
#include <TMVA/Reader.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
// analysis utilities
#include "TMVAClusterParameters.hxx"
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
  bool        do_read_cut;  // apply cuts while reading ntuple
}  DefaultOptions = {
  "./input/forNewTrainingMacro_noNonzeroEvts_andDefinitePrimary.evt5Ke210pim_central.d14m9y2024.root",
  "ntForCalib",
  "test.root",
  "tmva_test",
  "TMVARegression",
  true,
  false
};



// ============================================================================
//! Train and apply a TMVA model for BHCal cluster calibration
// ============================================================================
void TrainAndApplyBHCalClusterCalibration(const Options& opt = DefaultOptions) {

  // grab calculation parameters
  TMVAHelper::Parameters param = TMVAClusterParameters::GetParameters();

  // lower verbosity & announce start
  gErrorIgnoreLevel = kError;
  std::cout << "\n  Beginning calibration training and evaluation macro..." << std::endl;

  // --------------------------------------------------------------------------
  // Open input/outputs
  // --------------------------------------------------------------------------

  // open files
  TFile* inToTrain = new TFile(opt.in_file.data(),  "read");
  TFile* inToApply = new TFile(opt.in_file.data(),  "read");
  TFile* output    = new TFile(opt.out_file.data(), "recreate");
  if (!inToTrain || !inToApply || !output) {
    std::cerr << "PANIC: couldn't open a file!\n"
              << "       input  = " << inToTrain << "\n"
              << "       output = " << output
              << std::endl;
    assert(output && inToTrain && inToApply);
  }

  // print input/output files
  std::cout << "    Opened input/output files:\n"
            << "      input file  = " << opt.in_file << "\n"
            << "      output file = " << opt.out_file
            << std::endl;

  // grab input tuple
  TNtuple* ntToTrain = (TNtuple*) inToTrain -> Get(opt.in_tuple.data());
  TNtuple* ntToApply = (TNtuple*) inToApply -> Get(opt.in_tuple.data());
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

  // create tmva helpers
  TMVAHelper::Trainer train_helper( param.variables, param.methods );
  TMVAHelper::Reader  read_helper( param.variables, param.methods );
  train_helper.SetFactoryOptions(param.opts_factory);
  train_helper.SetTrainOptions(param.opts_training);
  read_helper.SetOptions(param.opts_reading);
  std::cout << "    Create TMVA helpers." << std::endl;

  // collect input leaves into a single vector
  std::vector<std::string> inputs;
  for (const auto& useAndVar : param.variables) {
    inputs.push_back(useAndVar.second);
  }

  // create input/output helpers
  NTupleHelper in_helper( inputs );
  NTupleHelper out_helper( read_helper.GetOutputs() );

  // set input/output tuple branches
  TNtuple* ntOutput = new TNtuple("ntTmvaOutput", "Output of TMVA regression", out_helper.CompressVariables().data());
  in_helper.SetBranches(ntToApply);
  std::cout << "    Set input/output tuple branches." << std::endl;


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
  train_helper.LoadVariables(loader, param.add_spectators);
  std::cout << "      Loaded variables..." << std::endl;

  // add tree & prepare for training
  loader -> AddRegressionTree(ntToTrain, param.tree_weight);
  loader -> PrepareTrainingAndTestTree(param.training_cuts, train_helper.CompressTrainingOptions().data());
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
  // Apply tmva models
  // --------------------------------------------------------------------------

  // instantiate formula object for applying ntuple cuts
  TTreeFormula* selector = new TTreeFormula("selector", param.reading_cuts, ntToApply);

  // instantiate reader
  TMVA::Reader* reader = new TMVA::Reader(read_helper.CompressOptions().data());
  std::cout << "    Begin applying calibration models:" << std::endl;

  // add input variables to reader, book methods
  read_helper.ReadVariables(reader, in_helper);
  read_helper.BookMethodsToRead(reader, opt.out_tmva, opt.name_tmva);
  std::cout << "      Added variables and methods to read." << std::endl;

  // get number of events for application
  const uint64_t nEntries = ntToApply -> GetEntries();
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
    const uint64_t bytes = ntToApply -> GetEntry(iEntry);
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

    // apply cuts if need be
    const bool isInCut = selector -> EvalInstance();
    if (opt.do_read_cut && !isInCut) continue;

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
  output    -> cd();
  ntOutput  -> Write(); 
  output    -> Close();
  inToTrain -> cd();
  inToTrain -> Close();
  inToApply -> cd();
  inToApply -> Close();

  // delete tmva objects
  delete factory;
  delete loader;
  delete reader;

  // announce end & exit
  std::cout << "  Finished BHCal calibration script!\n" << std::endl;
  return;

}

// end ========================================================================
