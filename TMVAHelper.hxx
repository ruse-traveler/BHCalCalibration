/// ===========================================================================
/*! \file   TMVAHelper.hxx
 *  \author Derek Anderson
 *  \date   09.11.2024
 *
 *  A lightweight class to help work with TMVA
 */
/// ===========================================================================

#ifndef TMVAHelper_hxx
#define TMVAHelper_hxx

// c++ utilities
#include <map>
#include <string>
#include <vector>
#include <stdio.h>
#include <cassert>
#include <utility>
#include <iostream>
// root libraries
#include <TString.h>
// tmva components
#include <TMVA/Tools.h>
#include <TMVA/Types.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
// analysis utilities
#include "NTupleHelper.h"



// ============================================================================
//! TMVA Helper
// ============================================================================
/*! A small class to help work with
 *  ROOT TMVA.
 */
class TMVAHelper {

  public:

    // ------------------------------------------------------------------------
    //! Use cases for input variables
    // ------------------------------------------------------------------------
    enum Use {Target, Train, Watch};

    // ------------------------------------------------------------------------
    //! Map of method names to algorithm type
    // ------------------------------------------------------------------------
    /*! This map enables automation of booking methods: some variants
     *  of algorithms (e.g. BDTG) cause TMVA to emit a std::runtime_error
     *  when attempting to look up the type based on the method name via
     *  the TMVA::Types helper object.
     */ 
    std::map<std::string, TMVA::Types::EMVA> MapNameToType = {
      {"BDT",      TMVA::Types::EMVA::kBDT},
      {"BDTG",     TMVA::Types::EMVA::kBDT},
      {"DNN",      TMVA::Types::EMVA::kDNN},
      {"FDA_GA",   TMVA::Types::EMVA::kFDA},
      {"FDA_GAMT", TMVA::Types::EMVA::kFDA},
      {"FDA_MC",   TMVA::Types::EMVA::kFDA},
      {"FDA_MT",   TMVA::Types::EMVA::kFDA},
      {"KNN",      TMVA::Types::EMVA::kKNN},
      {"LD",       TMVA::Types::EMVA::kLD},
      {"MLP",      TMVA::Types::EMVA::kMLP},
      {"PDEFoam",  TMVA::Types::EMVA::kPDEFoam},
      {"PDERS",    TMVA::Types::EMVA::kPDERS},
      {"SVM",      TMVA::Types::EMVA::kSVM}
    };

  private:

    // tmva option members
    std::vector<std::string> m_vecFactoryOpt;
    std::vector<std::string> m_vecTrainOpt;
    std::vector<std::string> m_vecReadOpt;

    // method members
    std::vector<std::string> m_methods;

    // input members
    std::vector<std::string> m_watchers;
    std::vector<std::string> m_trainers;
    std::vector<std::string> m_targets;

    // output members
    std::vector<std::string> m_outputs;

    // ------------------------------------------------------------------------
    //! Set input variables
    // ------------------------------------------------------------------------
    inline void SetInputVariables(const std::vector<std::pair<Use, std::string>>& inputs) {

      // assign input leaf to relevant vector
      for (const auto& input : inputs) {
        switch (input.first) {

          case Use::Target:
            m_targets.push_back( input.second );
            break;

          case Use::Train:
            m_trainers.push_back( input.second );
            break;

          case Use::Watch:
            [[fallthrough]];

          default:
            m_watchers.push_back( input.second );
            break;

        }
      }  // end input loop
      return;

    }  // end 'SetInputVariables(std::vector<std::pair<Use, std::string>>&)'

    // ------------------------------------------------------------------------
    //! Set methods to use
    // ------------------------------------------------------------------------
    //   - TODO might need a map and another vector for method-specific options
    inline void SetMethods(const std::vector<std::string>& methods) {

      m_methods = methods;
      return;

    }  // end 'SetMethods(std::vector<std::string>&)'

    // ------------------------------------------------------------------------
    //! Generate list of targets
    // ------------------------------------------------------------------------
    inline void GenerateRegressionOutputs() {

      for (const std::string& method : m_methods) {
        for (const std::string& target : m_targets) {
          m_outputs.push_back(target + "_" + method);
        }
      }
      return;

    }  // end 'GenerateListOfTargers()'

    // ------------------------------------------------------------------------
    //! Helper method to compress vector of strings into a colon-separated list
    // ------------------------------------------------------------------------
    inline std::string CompressList(const std::vector<std::string>& strings) const {

      std::string compressed("");
      for (std::size_t iString = 0; iString < strings.size(); ++iString) {
        compressed.append(strings[iString]);
        if (iString + 1 < strings.size()) {
          compressed.append(":");
        }
      }
      return compressed;

    }  // end 'CompressList(std::vector<std::string>&)'

  public:

    // ------------------------------------------------------------------------
    //! Getters
    // ------------------------------------------------------------------------
    inline std::vector<std::string> GetFactoryOptions()  const {return m_vecFactoryOpt;}
    inline std::vector<std::string> GetTrainingOptions() const {return m_vecTrainOpt;}
    inline std::vector<std::string> GetReadingOptions()  const {return m_vecReadOpt;}
    inline std::vector<std::string> GetMethods()         const {return m_methods;}
    inline std::vector<std::string> GetSpectatingVars()  const {return m_watchers;}
    inline std::vector<std::string> GetTrainingVars()    const {return m_trainers;}
    inline std::vector<std::string> GetTargetVars()      const {return m_targets;}
    inline std::vector<std::string> GetOutputs()         const {return m_outputs;}

    // ------------------------------------------------------------------------
    //! Compress lists into colon-separated ones
    // ------------------------------------------------------------------------
    inline std::string CompressReadOptions()    const {return CompressList(m_vecReadOpt);}
    inline std::string CompressTrainOptions()   const {return CompressList(m_vecTrainOpt);}
    inline std::string CompressFactoryOptions() const {return CompressList(m_vecFactoryOpt);}

    // ------------------------------------------------------------------------
    //! Setters
    // ------------------------------------------------------------------------
    inline void SetReadOptions(const std::vector<std::string>& options)    {m_vecReadOpt    = options;}
    inline void SetTrainOptions(const std::vector<std::string>& options)   {m_vecTrainOpt   = options;}
    inline void SetFactoryOptions(const std::vector<std::string>& options) {m_vecFactoryOpt = options;}

    // ------------------------------------------------------------------------
    //! Add variables to data loader
    // ------------------------------------------------------------------------
    inline void LoadVariables(TMVA::DataLoader* loader, const bool add_watchers = false) {

      // add targets
      for (const std::string& target : m_targets) {
        loader -> AddTarget(target);
      }

      // add training variables
      for (const std::string& train : m_trainers) {
        loader -> AddVariable(train);
      }

      // if needed, add spectators
      if (add_watchers) {
        for (const std::string& spec : m_watchers) {
          loader -> AddSpectator(spec);
        }
      }
      return;

    }  // end 'LoadVariables(TMVA::DataLoader*, bool)'

    // ------------------------------------------------------------------------
    //! Add NTuple variables to reader
    // ------------------------------------------------------------------------
    inline void ReadVariables(TMVA::Reader* reader, NTupleHelper& helper) {

      for (const std::string& train : m_trainers) {
        if (!helper.m_index.count(train)) {
          std::cerr << "WARNING: trying to add variable '" << train << "' which is not in input NTuple." << std::endl;
          continue;
        } else {
          reader -> AddVariable(train.data(), helper.m_values.at(m_index[train]));
        }
      }
      return;

    }  // end 'ReadVariables(TMVA::Reader*)'

    // ------------------------------------------------------------------------
    //! Book methods to train
    // ------------------------------------------------------------------------
    inline void BookMethodsToTrain(TMVA::Factory* factory, TMVA::DataLoader* loader) {

      // book each method currently set
      for (const std::string& method : m_methods) {
        factory -> BookMethod(
          loader,
          MapNameToType[method],
          method.data()
        );
      }
      return;

    }  // end 'BookMethodsToTrain(TMVA::Factory*, TMVA::DataLoader*)'

    // ------------------------------------------------------------------------
    //! Book methods to read by providing the path to a directory
    // ------------------------------------------------------------------------
    /*! For each method stored, helper will check for a weights file with
     *  a name conforming to "<method>.weights.xml" in the directory specified
     *  by `directory`. If the file exists, helper will book the method to be
     *  evaulated with that file.
     */
    inline void BookMethodsToRead(TMVA::Reader* reader, const std::string& directory) {

      // lambda to check if a file exists or not
      auto exists = [](const std::string& name) {
        if (FILE* file = fopen(name.c_str(), "r")) {
          fclose(file);
          return true;
        } else {
          return false;
        }
      };

      for (const std::string& method : m_methods) {

        // construct full path, skip if file does not exist
        const std::string path = directory + "/" method + ".weights.xml";
        if (!exists(path)) {
          std::cerr << "WARNING: file '" << path << "' doesn't exist! Not booking method!" << std::endl;
          continue;
        }

        // otherwise, construct title and book method
        const std::string title = method + " method";
        reader -> BookMVA(title, path);

      }
      return;

    }  // end 'BookMethodsToRead(TMVA::Reader*, directory)'

    // ------------------------------------------------------------------------
    //! Default ctor/dtor
    // ------------------------------------------------------------------------
    TMVAHelper()  {};
    ~TMVAHelper() {};

    // ------------------------------------------------------------------------
    //! ctor accepting a list of variable-use pairs and a list of methods
    // ------------------------------------------------------------------------
    TMVAHelper(
      const std::vector<std::pair<Use, std::string>>& inputs,
      const std::vector<std::string>& vecMethods
    ) {

      // set inputs and methods
      SetInputVariables(inputs);
      SetMethods(vecMethods);

      // then set the corresponding outputs
      GenerateRegressionOutputs();

    }  // end ctor(std::vector<std::pair<Use, std::string>>&, std::vector<std::string>&)'

};  // end TMVAHelper

#endif

// end ========================================================================
