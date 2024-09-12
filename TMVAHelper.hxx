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
#include <cassert>
#include <utility>
// root libraries
#include <TString.h>
// tmva components
#include <TMVA/Tools.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>



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

  private:

    // input members
    std::vector<std::string> m_watchers;
    std::vector<std::string> m_trainers;
    std::vector<std::string> m_targets;

    // output members
    std::map<std::string, std::size_t> m_index; 
    std::vector<float>                 m_values;

    // method members
    std::vector<std::string> m_methods;

    // tmva option members
    std::vector<std::string> m_vecFactoryOpt;
    std::vector<std::string> m_vecTrainOpt;
    std::vector<std::string> m_vecReadOpt;

    // ------------------------------------------------------------------------
    //! Set input variables
    // ------------------------------------------------------------------------
    inline void SetInputVariables(const std::vector<std::pair<Use, std::string>>& vecInput) {

      // assign input leaf to relevant vector
      for (const auto& input : vecInput) {
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
      //   - TODO might be useful to move map generation to another function
      std::size_t iTrain = 0;
      for (const std::string& train : m_trainers) {
        loader -> AddVariable(train);
        m_index[train] = iTrain;
        ++iTrain;
      }
      m_values.resize(m_trainers.size());

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

      // for looking up method type based on name
      TMVA::Types types = TMVA::Types::Instance();

      // book each method currently set
      for (const std::string& method : m_methods) {
        factory -> BookMethod(
          loader,
          types.GetMethodType( TString(method.data()) ),
          method.data()
        );
      }

      // destory types and exit
      types.DestroyInstance();
      return;

    }  // end 'BookMethodsToTrain(TMVA::Factory*, TMVA::DataLoader*)'

    // ------------------------------------------------------------------------
    //! Default ctor/dtor
    // ------------------------------------------------------------------------
    TMVAHelper()  {};
    ~TMVAHelper() {};

    // ------------------------------------------------------------------------
    //! ctor accepting a list of variable-use pairs and a list of methods
    // ------------------------------------------------------------------------
    TMVAHelper(
      const std::vector<std::pair<Use, std::string>>& vecInput,
      const std::vector<std::string>& vecMethods
    ) {

      SetInputVariables(vecInput);
      SetMethods(vecMethods);

      /* TODO need to generate vector of output targets */

    }  // end ctor(std::vector<std::pair<Use, std::string>>&, std::vector<std::string>&)'

};  // end TMVAHelper

#endif

// end ========================================================================
