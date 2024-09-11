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
#include <utility>
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

    // ------------------------------------------------------------------------
    //! Map of TMVA method name onto type enum
    // ------------------------------------------------------------------------
    const std::map<std::string, int> MapNameToIndex = {
      {"PDERS",    TMVA::Types::kPDERS},
      {"PDEFoam",  TMVA::Types::kPDEFoam},
      {"KNN",      TMVA::Types::kKNN},
      {"LD",       TMVA::Types::kLD},
      {"FDA_MC",   TMVA::Types::kFDA},
      {"FDA_GA",   TMVA::Types::kFDA},
      {"FDA_GAMT", TMVA::Types::kFDA},
      {"MLP",      TMVA::Types::kMLP},
      {"DNN_CPU",  TMVA::Types::kDL},
      {"DNN_GPU",  TMVA::Types::kDL},
      {"SVM",      TMVA::Types::kSVM},
      {"BDT",      TMVA::Types::kBDT},
      {"BDTG",     TMVA::Types::kBDT}
    };

  private:

    // input members
    std::vector<std::string> m_vecSpecLeaf;
    std::vector<std::string> m_vecTrainLeaf;
    std::vector<std::string> m_vecTargetLeaf;

    // tmva option members
    std::vector<std::string> m_vecFactoryOpt;
    std::vector<std::string> m_vecTrainerOpt;
    std::vector<std::string> m_vecReaderOpt;

  public:

    // ------------------------------------------------------------------------
    //! Set input variables
    // ------------------------------------------------------------------------
    inline void SetInputVariables(const std::vector<std::pair<Use, std::string>>& vecInput) {

      // assign input leaf to relevant vector
      for (const auto& input : vecInput) {
        switch (input.first) {

          case Use::Target:
            m_vecTargetLeaf.push_back( input.second );
            break;

          case Use::Train:
            m_vecTrainLeaf.push_back( input.second );
            break;

          case Use::Watch:
            [[fallthrough]];

          default:
            m_vecSpecLeaf.push_back( input.second );
            break;

        }
      }  // end input loop
      return;

    }  // end 'SetInputVariables(std::vector<std::pair<std::string, std::string>>&)'

    // ------------------------------------------------------------------------
    //! Default ctor/dtor
    // ------------------------------------------------------------------------
    TMVAHelper()  {};
    ~TMVAHelper() {};

};  // end TMVAHelper

#endif

// end ========================================================================
