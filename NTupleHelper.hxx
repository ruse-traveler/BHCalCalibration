/// ===========================================================================
/*! \file   NTupleHelper.hxx
 *  \author Derek Anderson
 *  \date   09.09.2024
 *
 *  A lightweight class to help work with TNTuple's.
 */
/// ===========================================================================

#ifndef NTupleHelper_hxx
#define NTupleHelper_hxx

// c++ utilities
#include <map>
#include <limits>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
// root libraries
#include <TNtuple.h>

// forward declaration of TMVAHelper::Reader
namespace TMVAHelper {
  class Reader;
}



// ============================================================================
//! NTuple Helper
// ============================================================================
/*! A small class to help work with
 *  ROOT TNtuple's.
 */
class NTupleHelper {

  private:

    // data members
    std::vector<float>                 m_values;
    std::vector<std::string>           m_variables;
    std::map<std::string, std::size_t> m_index;

  public:

    // ------------------------------------------------------------------------
    //! Getters
    // ------------------------------------------------------------------------
    inline std::vector<float>       GetValues()    const {return m_values;}
    inline std::vector<std::string> GetVariables() const {return m_variables;}

    // ------------------------------------------------------------------------
    //! Get a specific variable
    // ------------------------------------------------------------------------
    inline float GetVariable(const std::string& var) {

      // check if variable exists
      if (!m_index.count(var)) {
        assert(m_index.count(var));
      }

      // then get variable
      return m_values.at(m_index[var]);

    }  // end 'GetVariable(std::string&)'

    // ------------------------------------------------------------------------
    //! Set a variable
    // ------------------------------------------------------------------------
    inline void SetVariable(const std::string& var, const float val) {

      // check if variable exists
      if (!m_index.count(var)) {
        assert(m_index.count(var));
      }

      // then set variable
      m_values.at(m_index[var]) = val;
      return;

    }  // end 'SetVariable(std::string&, float)'

    // ------------------------------------------------------------------------
    //! Assign variables to TNtuple branches
    // ------------------------------------------------------------------------
    inline void SetBranches(TNtuple* tuple) {

      for (const std::string& var : m_variables) {
        tuple -> SetBranchAddress(var.data(), &m_values.at(m_index[var]));
      }
      return;

    }  // end 'SetBranches(TNtuple*)'

    // ------------------------------------------------------------------------
    //! Reset values
    // ------------------------------------------------------------------------
    inline void ResetValues() {

      std::fill(
        m_values.begin(),
        m_values.end(),
        -1. * std::numeric_limits<float>::max()
      );
      return;      

    }  // end 'ResetValues()'

    // ------------------------------------------------------------------------
    //! Compress list of variables into a colon-seperated list
    // ------------------------------------------------------------------------
    inline std::string CompressVariables() const {

      std::string compressed("");
      for (std::size_t iVar = 0; iVar < m_variables.size(); ++iVar) {
        compressed.append(m_variables[iVar]);
        if (iVar + 1 < m_variables.size()) {
          compressed.append(":");
        }
      }
      return compressed;

    }  // end 'CompressVariables()'

    // ------------------------------------------------------------------------
    //! Default ctor/dtor
    // ------------------------------------------------------------------------
    NTupleHelper()  {};
    ~NTupleHelper() {};

    // ------------------------------------------------------------------------
    //! ctor accepting a list of variables
    // ------------------------------------------------------------------------
    NTupleHelper(const std::vector<std::string>& vars) {

      // set variables
      m_variables = vars;

      // create map between name and position in vector
      std::size_t iVar = 0;
      for (const std::string& var : m_variables) {
        m_index[var] = iVar;
        ++iVar;
      }
      m_values.resize(m_variables.size());

    }  // end ctor(std::vector<std::string>&)'

    // make TMVAHelper::Reader a friend
    friend class TMVAHelper::Reader;

};  // end NTupleHelper

#endif

// end ========================================================================
