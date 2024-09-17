/// ===========================================================================
/*! \file   TMVAClusterParameters.hxx
 *  \author Derek Anderson
 *  \date   09.17.2024
 *
 *  A small namespace to define the parameters
 *  TMVA uses to train/apply various ML models
 *  for the calibration of energy for clusters
 *  in the BHCal and BIC.
 */
/// ===========================================================================

#ifndef TMVAClusterParameters_hxx
#define TMVAClusterParameters_hxx

// c++ utilities
#include <map>
#include <string>
#include <vector>
#include <utility>
// tmva components
#include <TMVA/Types.h>
// analysis utilities
#include "TMVAHelper.hxx"



// ============================================================================
//! Cluster TMVA Parameters
// ============================================================================
/*! TODO clean up a couple places to make slightly more readable
 */ 
namespace TMVAClusterParameters {

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

  // methods to use and options
  const std::vector<std::pair<std::string, std::string>> vecMethodAndOpt = {
    {
      "LD",
      "!H:!V:VarTransform=None"
    },
    {
      "KNN",
      "nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim"
    },
    {
      "MLP",
      "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator"
    },
    {
      "BDTG",
      "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:MaxDepth=4"
    },
    {
      "FDA_GA",
      "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:PopSize=100:Cycles=3:Steps=30:Trim=True:SaveBestGen=1:VarTransform=Norm"
    }
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
  const float treeWeight(1.0);
  const TCut  trainCut("(eSumBHCal>=0)&&(eSumBEMC>=0)&&(abs(hLeadBHCal)<1.1)&&(abs(hLeadBEMC)<1.1)");
  const TCut  readCut("(eLeadBEMC>0.5)&&(eLeadBEMC<100)");



  // --------------------------------------------------------------------------
  //! Collect options into parameter struct
  // --------------------------------------------------------------------------
  TMVAHelper::Parameters GetParameters() {

    TMVAHelper::Parameters param {
      .variables      = vecUseAndVar,
      .methods        = vecMethodAndOpt,
      .opts_factory   = vecFactoryOpts,
      .opts_training  = vecTrainOpts,
      .opts_reading   = vecReadOpts,
      .add_spectators = addSpectators,
      .training_cuts  = trainCut,
      .reading_cuts   = readCut,
      .tree_weight    = treeWeight
    };
    return param;

  }  // end 'GetParameters()'

}  // end TMVAClusterParameters namespace

#endif

// end ========================================================================
