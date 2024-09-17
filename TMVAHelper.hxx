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
#include <algorithm>
// root libraries
#include <TString.h>
// tmva components
#include <TMVA/Tools.h>
#include <TMVA/Types.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
// analysis utilities
#include "NTupleHelper.hxx"



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
    std::vector<bool>        m_read;
    std::vector<std::string> m_methods;

    // input members
    std::vector<std::string> m_watchers;
    std::vector<std::string> m_trainers;
    std::vector<std::string> m_targets;

    // output members
    std::vector<float>                 m_out_values;
    std::vector<std::string>           m_out_variables;
    std::map<std::string, std::size_t> m_out_index;

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
      m_read.resize( m_methods.size(), true );
      return;

    }  // end 'SetMethods(std::vector<std::string>&)'

    // ------------------------------------------------------------------------
    //! Generate list of regression outputs
    // ------------------------------------------------------------------------
    /*! Will generate list of targets and the regression
     *  outputs. Every specified target will be evaluated
     *  for every specified method.
     */
    inline void GenerateRegressionOutputs() {

      // 1st load targets
      std::size_t iOut = 0;
      for (const std::string& target : m_targets) {
        m_out_variables.push_back(target);
        m_out_index[m_out_variables.back()] = iOut;
        ++iOut;
      }

      // then generate list of regression outputs
      for (const std::string& method : m_methods) {
        for (const std::string& target : m_targets) {
          m_out_variables.push_back(target + "_" + method);
          m_out_index[m_out_variables.back()] = iOut;
          ++iOut;
        }
      }

      // resize values vector acordingly & exit
      m_out_values.resize(iOut);
      return;

    }  // end 'GenerateListOfTargers()'

    // ------------------------------------------------------------------------
    //! Helper method to check if a file exists or not
    // ------------------------------------------------------------------------
    inline bool DoesFileExist(const std::string& name) const {

      bool exists;
      if (FILE* file = fopen(name.c_str(), "r")) {
        fclose(file);
        exists = true;
      } else {
        exists = false;
      }
      return exists;

    }  // end 'DoesFileExist(std::string)'

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
    //! Compress lists into colon-separated ones
    // ------------------------------------------------------------------------
    inline std::string CompressReadOptions()    const {return CompressList(m_vecReadOpt);}
    inline std::string CompressTrainOptions()   const {return CompressList(m_vecTrainOpt);}
    inline std::string CompressFactoryOptions() const {return CompressList(m_vecFactoryOpt);}

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
    inline std::vector<std::string> GetOutVars()         const {return m_out_variables;}

    // ------------------------------------------------------------------------
    //! Get a specific output variable
    // ------------------------------------------------------------------------
    inline float GetVariable(const std::string& var) {

      // check if variable exists
      if (!m_out_index.count(var)) {
        assert(m_out_index.count(var));
      }

      // then get variable
      return m_out_values.at(m_out_index[var]);

    }  // end 'GetVariable(std::string&)'

    // ------------------------------------------------------------------------
    //! Reset output values
    // ------------------------------------------------------------------------
    inline void ResetValues() {

      std::fill(
        m_out_values.begin(),
        m_out_values.end(),
        -1. * std::numeric_limits<float>::max()
      );
      return;

    }  // end 'ResetValues()'

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
    inline void ReadVariables(TMVA::Reader* reader, NTupleHelper& helper) const {

      for (const std::string& train : m_trainers) {
        if (!helper.m_index.count(train)) {
          std::cerr << "WARNING: trying to add variable '" << train << "' which is not in input NTuple." << std::endl;
          continue;
        } else {
          reader -> AddVariable(train.data(), &helper.m_values.at(helper.m_index[train]));
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
    inline void BookMethodsToRead(TMVA::Reader* reader, const std::string& directory, const std::string& name) {

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

    // ------------------------------------------------------------------------
    //! Book methods to read by providing a list of weight files
    // ------------------------------------------------------------------------
    /*! For file passed, helper wil check if file exists and, if it does,
     *  book the corresponding method to be evaluated with that file. This
     *  assumes that the insertion order is the same between the methods
     *  list and the provided list of files.
     */
    inline void BookMethodsToRead(TMVA::Reader* reader, const std::vector<std::string>& files) {

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

    // ------------------------------------------------------------------------
    //! Evaluate all booked methods
    // ------------------------------------------------------------------------
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
          m_out_values.at(m_out_index[output]) = targets.at(iTarget);
        }
      }

      // then collect training targets in output
      for (const std::string& target : m_targets) {
        m_out_values.at(m_out_index[target]) = helper.GetVariable( target );
      }
      return;

    }  // end 'EvaluateMethods(TMVA::Reader*, NTupleHelper&)'

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
