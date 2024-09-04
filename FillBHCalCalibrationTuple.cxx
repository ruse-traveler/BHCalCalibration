/// ===========================================================================
/*! \file   FillBHCalCalibrationTuple.cxx
 *  \author Derek Anderson
 *  \date   09.04.2024
 *
 *  A ROOT macro to read EICrecon output (either `*.podio.root` or
 *  `*.tree.edm4eic.root` and fill an NTuple for training a ML model.
 */
/// ===========================================================================

#define FillBHCalCalibrationTuple_cxx

// c++ utilities
#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <optional>
// root libraries
#include <TFile.h>
#include <TNtuple.h>
#include <TSystem.h>
// podio libraries
#include <podio/Frame.h>
#include <podio/CollectionBase.h>
#include <podio/ROOTFrameReader.h>
// edm4eic types
#include <edm4eic/ClusterCollection.h>
#include <edm4eic/ReconstructedParticleCollection.h>

  // TODO REMOVE ME
  // global constants
  enum CONST {
    NCalibVars  = 51,
    NSciFiLayer = 12,
    NImageLayer = 4,
    NRange      = 2,
    NComp       = 3
  };



// ============================================================================
//! Struct to consolidate user options
// ============================================================================
struct Options {
  std::string in_file;      // input file
  std::string out_file;     // output file
  std::string gen_par;      // generated particles
  std::string hcal_clust;   // hcal cluster collection
  std::string ecal_clust;   // ecal (scfi + imaging) cluster collection
  std::string scfi_clust;   // ecal (scfi) cluster collection
  std::string image_clust;  // ecal (imaging) cluster/layer collection
  bool        do_progress;  // print progress through frame loop
} DefaultOptions = {
  "./input/testNewSplitMergeEdits_useSetForConsumed.e10pim_central.d2m9y2024.podio.root",
  "output.root",
  "GeneratedParticles",
  "HcalBarrelClusters",
  "EcalBarrelClusters",
  "EcalBarrelScFiClusters",
  "EcalBarrelImagingLayers",
  true
};



// ============================================================================
//! Small class to manage output
// ============================================================================
class Output {

  private:

    // data members
    std::vector<std::string> m_variables;
    std::vector<float>       m_values;

  public:

    // --------------------------------------------------------------------------
    //! Compress list of variables into a colon-seperated list
    // --------------------------------------------------------------------------
    inline std::string CompressVariables() const {

      std::string compressed("");
      for (
        std::size_t iVariable = 0;
        iVariable < m_variables.size();
        ++iVariable
      ) {
        compressed.append(m_variables[iVariable] );
        if (iVariable + 1 != m_variables.size()) {
          compressed.append(":");
        }
      }
      return compressed;

    }  // end 'CompressVariables()'

    // ------------------------------------------------------------------------
    // Access values
    // -----------------------------------------------------------------------
    inline std::vector<float> GetValues() const {

      return m_values;

    }  // end 'GetValues()'

    // ctor/dtor
    Output()  {};
    ~Output() {};

};  // end Output



// ============================================================================
//! Fill BHCal calibration NTuple
// ============================================================================
void FillBHCalCalibrationTuple(const Options& opt = DefaultOptions) {

  // announce start of macro
  std::cout << "\n  Beginning calibration tuple-filling macro!" << std::endl;

  // --------------------------------------------------------------------------
  // Open input/outputs
  // --------------------------------------------------------------------------

  // open file w/ frame reader
  podio::ROOTFrameReader reader = podio::ROOTFrameReader();
  reader.openFile( opt.in_file );

  // open output file
  TFile* output = new TFile(opt.out_file.data(), "recreate");
  if (!output) {
    std::cerr << "PANIC: couldn't open output file!" << std::endl;
    assert(output);
  }

  // print input/output
  std::cout << "    Opened input/output files:\n"
            << "      input file  = " << opt.in_file << "\n"
            << "      output file = " << opt.out_file
            << std::endl;

  // --------------------------------------------------------------------------
  // Loop over input frames
  // --------------------------------------------------------------------------

  // get no. of frames and announce start of frame loop
  const uint64_t nFrames = reader.getEntries(podio::Category::Event);
  std::cout << "    Starting frame loop: " << reader.getEntries(podio::Category::Event) << " frames to process." << std::endl;

  // iterate through frames
  for (uint64_t iFrame = 0; iFrame < nFrames; ++iFrame) {

    // announce progress
    if (opt.do_progress) {
      std::cout << "      Processing frame " << iFrame + 1 << "/" << nFrames << "...";
      if (iFrame + 1 < nFrames) {
        std::cout << "\r" << std::flush;
      } else {
        std::cout << std::endl;
      }
    }

    // grab frame
    auto frame = podio::Frame( reader.readNextEntry(podio::Category::Event) );

    // grab needed collections
    auto& genParticles  = frame.get<edm4eic::ReconstructedParticleCollection>( opt.gen_par );
    auto& hcalClusters  = frame.get<edm4eic::ClusterCollection>( opt.hcal_clust );
    auto& ecalClusters  = frame.get<edm4eic::ClusterCollection>( opt.ecal_clust );
    auto& scfiClusters  = frame.get<edm4eic::ClusterCollection>( opt.scfi_clust );
    auto& imageClusters = frame.get<edm4eic::ClusterCollection>( opt.image_clust );

    // ------------------------------------------------------------------------
    // particle loop
    // ------------------------------------------------------------------------

     for (edm4eic::ReconstructedParticle particle : genParticles) {

       /* TODO fill in */

     }  // end particle loop 

    // ------------------------------------------------------------------------
    // hcal cluster loop
    // ------------------------------------------------------------------------

    for (edm4eic::Cluster hClust : hcalClusters) {

      /* TODO fill in */

    }  // end hcal cluster loop

    // ------------------------------------------------------------------------
    // ecal cluster loops
    // ------------------------------------------------------------------------

    // loop over combined ecal clusters
    for (edm4eic::Cluster eClust : ecalClusters) {

      /* TODO fill in */

    }  // end combined ecal cluster loop

    // loop over scfi ecal clusters
    for (edm4eic::Cluster sClust : scfiClusters) {

      /* TODO fill in */

    }  // end scfi cluster loop

    // loop over imaging ecal clusters (layers)
    for (edm4eic::Cluster iClust : imageClusters) {

    }  // end imaging cluster loop

    // ------------------------------------------------------------------------
    // fill ntuple
    //  -----------------------------------------------------------------------

    /* TODO goes here */

  }  // end frame loop
  std::cout << "    Finished frame loop" << std::endl;

  // save output & close files
  output -> cd();
  output -> Close();

  // announce end & exit
  std::cout << "  End of macro!\n" << std::endl;
  return;

}

// end ========================================================================
