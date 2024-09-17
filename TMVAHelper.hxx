/// ===========================================================================
/*! \file   TMVAHelper.hxx
 *  \author Derek Anderson
 *  \date   09.11.2024
 *
 *  A lightweight namespace to help work with TMVA
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
#include <algorithm>
// root libraries
#include <TCut.h>
#include <TString.h>
// tmva components
#include <TMVA/Tools.h>
#include <TMVA/Types.h>
#include <TMVA/Reader.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
// analysis utilities
#include "NTupleHelper.hxx"



// ============================================================================
//! TMVA Helper
// ============================================================================
/*! A small namespace to help work with
 *  ROOT TMVA.
 */
namespace TMVAHelper {

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
  inline std::map<std::string, TMVA::Types::EMVA> MapNameToType() {

    static std::map<std::string, TMVA::Types::EMVA> mapNameToType = {
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
    return mapNameToType;

  }  // end 'MapNameToType()'



  // --------------------------------------------------------------------------
  //! Helper method to check if a file exists or not
  // --------------------------------------------------------------------------
  inline bool DoesFileExist(const std::string& name) {

    if (FILE* file = fopen(name.c_str(), "r")) {
      fclose(file);
      return true;
    } else {
      return false;
    }

  }  // end 'DoesFileExist(std::string)'



  // --------------------------------------------------------------------------
  //! Helper method to compress vector of strings into a colon-separated list
  // --------------------------------------------------------------------------
  inline std::string CompressList(const std::vector<std::string>& strings) {

    std::string compressed("");
    for (std::size_t iString = 0; iString < strings.size(); ++iString) {
      compressed.append(strings[iString]);
      if (iString + 1 < strings.size()) {
        compressed.append(":");
      }
    }
    return compressed;

  }  // end 'CompressList(std::vector<std::string>&)'



  // ==========================================================================
  //! TMVA Parameters
  // ==========================================================================
  /*! Struct to consolidate training/application
   *  parameters used.
   */
  struct Parameters {

    // variables & methods
    std::vector<std::pair<Use, std::string>>         variables;  // input variables & usage
    std::vector<std::pair<std::string, std::string>> methods;    // methods to use & options

    // tmva options
    std::vector<std::string> opts_factory;   // general options
    std::vector<std::string> opts_training;  // training/testing options
    std::vector<std::string> opts_reading;   // reading/evaluating options

    // other parameters
    bool  add_spectators;  // whether or not to add spectator variables
    TCut  training_cuts;   // cuts to apply during training
    TCut  reading_cuts;    // cuts to apply during reading
    float tree_weight;     // weight of tree being trained on

  };  // end TMVAHelper::Parameters



  // ==========================================================================
  //! Base Helper
  // ==========================================================================
  /*! Base TMVA Helper class that the training and
   *  reading helpers inherit from. Defines common
   *  methods used betwee the two.
   */
  class Base {

    protected:

      // data members
      std::vector<std::string>           m_watchers;
      std::vector<std::string>           m_trainers;
      std::vector<std::string>           m_targets;
      std::vector<std::string>           m_methods;
      std::map<std::string, std::string> m_opts_method;

      // ----------------------------------------------------------------------
      //! Set methods to use
      // ----------------------------------------------------------------------
      inline void SetMethods(const std::vector<std::pair<std::string, std::string>>& methods) {

        for (const auto& methodAndOpt : methods) {
          m_methods.push_back( methodAndOpt.first );
          m_opts_method[methodAndOpt.first] = methodAndOpt.second;
        }
        return;

      }  // end 'SetMethods(std::vector<std::string>&)'

      // ----------------------------------------------------------------------
      //! Set input variables
      // ----------------------------------------------------------------------
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

    public:

      // ----------------------------------------------------------------------
      //! Getters
      // ----------------------------------------------------------------------
      inline std::vector<std::string> GetMethods()    const {return m_methods;}
      inline std::vector<std::string> GetSpectators() const {return m_watchers;}
      inline std::vector<std::string> GetTrainers()   const {return m_trainers;}
      inline std::vector<std::string> GetTargets()    const {return m_targets;}

      // ----------------------------------------------------------------------
      //! Default ctor/dtor
      // ----------------------------------------------------------------------
      Base()  {};
      ~Base() {};

  };  // end TMVAHelper::Base



  // ==========================================================================
  //! Training Helper
  // ==========================================================================
  /*! A small class to help training via
   *  models ROOT TMVA.
   */
  class Trainer : public Base {

    private:

      // data members
      std::vector<std::string> m_opts_factory;
      std::vector<std::string> m_opts_train;

    public:

      // ----------------------------------------------------------------------
      //! Compress lists into colon-separated ones
      // ----------------------------------------------------------------------
      inline std::string CompressFactoryOptions()  const {return CompressList(m_opts_factory);}
      inline std::string CompressTrainingOptions() const {return CompressList(m_opts_train);}

      // ----------------------------------------------------------------------
      //! Getters
      // ----------------------------------------------------------------------
      inline std::vector<std::string> GetFactoryOptions()  const {return m_opts_factory;}
      inline std::vector<std::string> GetTrainingOptions() const {return m_opts_train;}

      // ----------------------------------------------------------------------
      //! Setters
      // ----------------------------------------------------------------------
      inline void SetTrainOptions(const std::vector<std::string>& options)   {m_opts_train   = options;}
      inline void SetFactoryOptions(const std::vector<std::string>& options) {m_opts_factory = options;}

      // ----------------------------------------------------------------------
      //! Add variables to data loader
      // ----------------------------------------------------------------------
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
      //! Book methods to train
      // ------------------------------------------------------------------------
      inline void BookMethodsToTrain(TMVA::Factory* factory, TMVA::DataLoader* loader) {

        // book each method currently set
        for (const std::string& method : m_methods) {
          factory -> BookMethod(
            loader,
            MapNameToType()[method],
            method.data(),
            m_opts_method[method].data()
          );
        }
        return;

      }  // end 'BookMethodsToTrain(TMVA::Factory*, TMVA::DataLoader*)'

      // ----------------------------------------------------------------------
      //! Default ctor/dtor
      // ----------------------------------------------------------------------
      Trainer()  {};
      ~Trainer() {};

      // ----------------------------------------------------------------------
      //! ctor accepting a list of variable-use pairs and a list of methods
      // ----------------------------------------------------------------------
      Trainer(
        const std::vector<std::pair<Use, std::string>>& inputs,
        const std::vector<std::pair<std::string, std::string>>& methods
      ) {

        SetInputVariables(inputs);
        SetMethods(methods);

      }  // end ctor(std::vector<std::pair<Use, std::string>>&, std::vector<std::string>&)'

  };  // end TMVAHelper::Trainer



  // ==========================================================================
  //! Reading Helper
  // ==========================================================================
  /*! A small class to help reading/evaluating
   *  modles via ROOT TMVA.
   */
  class Reader : public Base {

    private:

      // data members
      std::vector<bool>                  m_read;
      std::vector<float>                 m_outvals;
      std::vector<std::string>           m_outvars;
      std::vector<std::string>           m_options;
      std::map<std::string, std::size_t> m_outdex;

      // ----------------------------------------------------------------------
      //! Generate list of regression outputs
      // ----------------------------------------------------------------------
      /*! Will generate list of targets and the regression
       *  outputs. Every specified target will be evaluated
       *  for every specified method.
       */
      inline void GenerateRegressionOutputs() {

        // first load targets
        std::size_t iOut = 0;
        for (const std::string& target : m_targets) {
          m_outvars.push_back( target );
          m_outdex[ m_outvars.back() ] = iOut;
          ++iOut;
        }

        // then generate list of regression outputs
        for (const std::string& method : m_methods) {
          for (const std::string& target : m_targets) {
            m_outvars.push_back( target + "_" + method );
            m_outdex[ m_outvars.back() ] = iOut;
            ++iOut;
          }
        }

        // resize values vector acordingly & exit
        m_outvals.resize(iOut);
        return;

      }  // end 'GenerateRegressionOutputs()'

    public:

      // ----------------------------------------------------------------------
      //! Compress list into colon-separated one
      // ----------------------------------------------------------------------
      inline std::string CompressOptions() const {return CompressList(m_options);}

      // ------------------------------------------------------------------------
      //! Setters
      // ------------------------------------------------------------------------
      inline void SetOptions(const std::vector<std::string>& options) {m_options = options;}

      // ----------------------------------------------------------------------
      //! Getters
      // ----------------------------------------------------------------------
      inline std::vector<std::string> GetOptions() const {return m_options;}
      inline std::vector<std::string> GetOutputs() const {return m_outvars;}

      // ----------------------------------------------------------------------
      //! Get a specific output variable
      // ----------------------------------------------------------------------
      inline float GetVariable(const std::string& var) {

        // check if variable exists
        if (!m_outdex.count(var)) {
          assert(m_outdex.count(var));
        }

        // then get variable
        return m_outvals.at(m_outdex[var]);

      }  // end 'GetVariable(std::string&)'

      // ----------------------------------------------------------------------
      //! Reset output values
      // ----------------------------------------------------------------------
      inline void ResetValues() {

        std::fill(
          m_outvals.begin(),
          m_outvals.end(),
          -1. * std::numeric_limits<float>::max()
        );
        return;

      }  // end 'ResetValues()'

      // ----------------------------------------------------------------------
      //! Add NTuple variables to reader
      // ----------------------------------------------------------------------
      inline void ReadVariables(TMVA::Reader* reader, NTupleHelper& helper) const {

        for (const std::string& train : m_trainers) {
          if (!helper.m_index.count(train)) {
            std::cerr << "WARNING: trying to add variable '" << train << "' which is not in input NTuple!" << std::endl;
            continue;
          } else {
            reader -> AddVariable(train.data(), &helper.m_values.at(helper.m_index[train]));
          }
        }
        return;

      }  // end 'ReadVariables(TMVA::Reader*, NTupleHelper&)'

      // ----------------------------------------------------------------------
      //! Book methods to read by providing the path to a directory
      // ----------------------------------------------------------------------
      /*! For each method stored, helper will check for a weights file with
       *  a name conforming to "<method>.weights.xml" in the directory specified
       *  by `directory`. If the file exists, helper will book the method to be
       *  evaulated with that file.
       */
      inline void BookMethodsToRead(
        TMVA::Reader* reader,
        const std::string& directory,
        const std::string& name
      ) {

        // reserve space for each method
        m_read.resize( m_methods.size(), true );

        // loop over all methods
        for (std::size_t iMethod = 0; iMethod < m_methods.size(); ++iMethod) {

          // construct full path
          const std::string path = directory + "/weights/" + name + "_" + m_methods[iMethod] + ".weights.xml";

          // skip if file does not exist
          if (!DoesFileExist(path)) {
            std::cerr << "WARNING: file '" << path << "' doesn't exist! Not booking method!" << std::endl;
            m_read.at(iMethod) = false;
            continue;
          }

          // otherwise, construct title and book method
          const std::string title = m_methods[iMethod] + " method";
          reader -> BookMVA(title, path);

        }  // end method loop
        return;

      }  // end 'BookMethodsToRead(TMVA::Reader*, std::string&, std::string&)'

      // ----------------------------------------------------------------------
      //! Book methods to read by providing a list of weight files
      // ----------------------------------------------------------------------
      /*! For file passed, helper wil check if file exists and, if it does,
       *  book the corresponding method to be evaluated with that file. This
       *  assumes that the insertion order is the same between the methods
       *  list and the provided list of files.
       */
      inline void BookMethodsToRead(
        TMVA::Reader* reader,
        const std::vector<std::string>& files
      ) {

        // reserve space for each method
        m_read.resize( m_methods.size(), true );

        // make sure input list has same dimension as method list
        if (files.size() != m_methods.size()) {
          assert(files.size() == m_methods.size());
        }

        // loop over provided files
        for (std::size_t iFile = 0; iFile < files.size(); ++iFile) {

          // skip if file does not exist
          if (!DoesFileExist(files[iFile])) {
            std::cerr << "WARNING: file '" << files[iFile] << "' doesn't exist! Not booking method '" << m_methods.at(iFile) <<"'!" << std::endl;
            m_read.at(iFile) = false;
            continue;
          }

          // otherwise, construct title and book method
          const std::string title = m_methods.at(iFile) + " method";
          reader -> BookMVA(title, files[iFile]);

        }  // end file loop
        return;

      }  // end 'BookMethodsToRead(TMVA::Reader*, std::vector<std::string>&)'

      // ----------------------------------------------------------------------
      //! Evaluate all booked methods
      // ----------------------------------------------------------------------
      inline void EvaluateMethods(TMVA::Reader* reader, NTupleHelper& helper) {

        // loop over all methods
        for (std::size_t iMethod = 0; iMethod < m_methods.size(); ++iMethod) {

          // if not evaluating method, continue
          if (!m_read[iMethod]) {
            continue;
          }

          // construct title & run evaluation
          const std::string        title   = m_methods[iMethod] + " method";
          const std::vector<float> targets = reader -> EvaluateRegression(title);

          // collect regression output
          for (std::size_t iTarget = 0; iTarget < m_targets.size(); ++iTarget) {
            const std::string output = m_targets[iTarget] + "_" + m_methods[iMethod];
            m_outvals.at(m_outdex[output]) = targets.at(iTarget);
          }
        }  // end method loop

        // then collect training targets in output
        for (const std::string& target : m_targets) {
          m_outvals.at(m_outdex[target]) = helper.GetVariable( target );
        }
        return;

      }  // end 'EvaluateMethods(TMVA::Reader*, NTupleHelper&)'

      // ----------------------------------------------------------------------
      //! Default ctor/dtor
      // ----------------------------------------------------------------------
      Reader()  {};
      ~Reader() {};

      // ----------------------------------------------------------------------
      //! ctor accepting a list of variable-use pairs and a list of methods
      // ----------------------------------------------------------------------
      Reader(
        const std::vector<std::pair<Use, std::string>>& inputs,
        const std::vector<std::pair<std::string, std::string>>& methods
      ) {

        SetInputVariables(inputs);
        SetMethods(methods);
        GenerateRegressionOutputs();

      }  // end ctor(std::vector<std::pair<Use, std::string>>&, std::vector<std::string>&)'

  };  // end TMVAHelper::Reader

}  // end TMVAHelper namespace

#endif

// end ========================================================================
