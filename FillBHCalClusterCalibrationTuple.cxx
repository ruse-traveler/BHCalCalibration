/// ===========================================================================
/*! \file   FillBHCalClusterCalibrationTuple.cxx
 *  \author Derek Anderson
 *  \date   09.04.2024
 *
 *  A ROOT macro to read EICrecon output (either `*.podio.root` or
 *  `*.tree.edm4eic.root` and fill an NTuple for training a ML model.
 */
/// ===========================================================================

#define FillBHCalClusterCalibrationTuple_cxx

// c++ utilities
#include <limits>
#include <string>
#include <vector>
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
#include <edm4eic/CalorimeterHitCollection.h>
#include <edm4eic/ReconstructedParticleCollection.h>
// edm4hep types
#include <edm4hep/Vector3f.h>
#include <edm4hep/utils/vector_utils.h>
// analysis utilities
#include "../../utility/NTupleHelper.hxx"



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
  std::string scfi_hits;    // ecal (scfi) hit collection
  std::string image_clust;  // ecal (imaging) cluster/layer collection
  std::string image_hits;   // ecal (imaging) hit collection
  bool        do_progress;  // print progress through frame loop
} DefaultOptions = {
  "./forNewCalibWorkflow.evt5Ke10pim_central.d14m9y2024.podio.root",
  "forNewTrainingMacro_noNonzeroEvts_andDefinitePrimary.evt5Ke10pim_central.d14m9y2024.root",
  "GeneratedParticles",
  "HcalBarrelClusters",
  "EcalBarrelClusters",
  "EcalBarrelScFiClusters",
  "EcalBarrelScFiRecHits",
  "EcalBarrelImagingLayers",
  "EcalBarrelImagingRecHits",
  true
};



// ============================================================================
//! Fill BHCal cluster calibration NTuple
// ============================================================================
void FillBHCalClusterCalibrationTuple(const Options& opt = DefaultOptions) {

  // --------------------------------------------------------------------------
  // calculation parameters
  // --------------------------------------------------------------------------

  // output variables
  NTupleHelper helper({
    "ePar",
    "fracParVsLeadBHCal",
    "fracParVsLeadBEMC",
    "fracParVsSumBHCal",
    "fracParVsSumBEMC",
    "fracLeadBHCalVsBEMC",
    "fracSumBHCalVsBEMC",
    "eLeadBHCal",
    "eLeadBEMC",
    "eSumBHCal",
    "eSumBEMC",
    "diffLeadBHCal",
    "diffLeadBEMC",
    "diffSumBHCal",
    "diffSumBEMC",
    "nHitsLeadBHCal",
    "nHitsLeadBEMC",
    "nClustBHCal",
    "nClustBEMC",
    "hLeadBHCal",
    "hLeadBEMC",
    "fLeadBHCal",
    "fLeadBEMC",
    "eLeadImage",
    "eSumImage",
    "eLeadScFi",
    "eSumScFi",
    "nClustImage",
    "nClustScFi",
    "hLeadImage",
    "hLeadScFi",
    "fLeadImage",
    "fLeadScFi",
    "eSumScFiLayer1",
    "eSumScFiLayer2",
    "eSumScFiLayer3",
    "eSumScFiLayer4",
    "eSumScFiLayer5",
    "eSumScFiLayer6",
    "eSumScFiLayer7",
    "eSumScFiLayer8",
    "eSumScFiLayer9",
    "eSumScFiLayer10",
    "eSumScFiLayer11",
    "eSumScFiLayer12",
    "eSumImageLayer1",
    "eSumImageLayer2",
    "eSumImageLayer3",
    "eSumImageLayer4",
    "eSumImageLayer5",
    "eSumImageLayer6"
  });

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

  // create output ntuple
  TNtuple* ntForCalib = new TNtuple("ntForCalib", "NTuple for calibration", helper.CompressVariables().c_str());

  // --------------------------------------------------------------------------
  // Loop over input frames
  // --------------------------------------------------------------------------
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
    auto& scfiHits      = frame.get<edm4eic::CalorimeterHitCollection>( opt.scfi_hits );
    auto& imageClusters = frame.get<edm4eic::ClusterCollection>( opt.image_clust );
    auto& imageHits     = frame.get<edm4eic::CalorimeterHitCollection>( opt.image_hits );

    // reset output values
    helper.ResetValues();

    // ------------------------------------------------------------------------
    // particle loop
    // ------------------------------------------------------------------------
    std::optional<edm4eic::ReconstructedParticle> optPrimary = std::nullopt;
    for (edm4eic::ReconstructedParticle particle : genParticles) {
      if (particle.getType() == 1) {
        optPrimary = particle;
        break;
      }
    }  // end particle loop

    // skip event if no primary found
    if (!optPrimary.has_value()) {
      continue;
    }
    edm4eic::ReconstructedParticle primary = optPrimary.value();

    // set particle output variables
    helper.SetVariable( "ePar", primary.getEnergy() );

    // ------------------------------------------------------------------------
    // hcal cluster loop
    // ------------------------------------------------------------------------
    edm4eic::Cluster hLeadClust;

    // find leading cluster, sum energies
    float eSumHCal  = 0.;
    float eLeadHCal = 0.;
    for (edm4eic::Cluster hClust : hcalClusters) {

      if (hClust.getEnergy() > eLeadHCal) {
        hLeadClust = hClust;
        eLeadHCal  = hClust.getEnergy();
      }
      eSumHCal += hClust.getEnergy();

    }  // end hcal cluster loop

    // fill lead hcal cluster variables
    helper.SetVariable( "eLeadBHCal", hLeadClust.getEnergy() );
    helper.SetVariable( "nHitsLeadBHCal", (float) hLeadClust.getHits().size() );
    helper.SetVariable( "hLeadBHCal", edm4hep::utils::eta(hLeadClust.getPosition()) );
    helper.SetVariable( "fLeadBHCal", edm4hep::utils::angleAzimuthal(hLeadClust.getPosition()) );

    // fill event-level output variables
    helper.SetVariable( "eSumBHCal", eSumHCal);
    helper.SetVariable( "nClustBHCal", (float) hcalClusters.size());
    helper.SetVariable( "fracParVsSumBHCal", eSumHCal / primary.getEnergy());
    helper.SetVariable( "fracParVsLeadBHCal", hLeadClust.getEnergy() / primary.getEnergy());
    helper.SetVariable( "diffSumBHCal", (eSumHCal - primary.getEnergy()) / primary.getEnergy());
    helper.SetVariable( "diffLeadBHCal", (hLeadClust.getEnergy() - primary.getEnergy()) / primary.getEnergy());

    // ------------------------------------------------------------------------
    // ecal (scfi + imaging) cluster loop
    // ------------------------------------------------------------------------
    edm4eic::Cluster eLeadClust;

    // loop over combined ecal clusters
    float eSumECal  = 0.;
    float eLeadECal = 0.;
    for (edm4eic::Cluster eClust : ecalClusters) {

      if (eClust.getEnergy() > eLeadECal) {
        eLeadClust = eClust;
        eLeadECal  = eClust.getEnergy();
      }
      eSumECal += eClust.getEnergy();

    }  // end combined ecal cluster loop

    // fill lead ecal cluster variables
    helper.SetVariable( "eLeadBEMC", eLeadClust.getEnergy() );
    helper.SetVariable( "nHitsLeadBEMC", (float) eLeadClust.getHits().size() );
    helper.SetVariable( "hLeadBEMC", edm4hep::utils::eta(eLeadClust.getPosition()) );
    helper.SetVariable( "fLeadBEMC", edm4hep::utils::angleAzimuthal(eLeadClust.getPosition()) );

    // fill event-level output variables
    helper.SetVariable( "eSumBEMC", eSumECal );
    helper.SetVariable( "nClustBEMC", (float) ecalClusters.size() );
    helper.SetVariable( "fracParVsSumBEMC", eSumECal / primary.getEnergy() );
    helper.SetVariable( "fracParVsLeadBEMC", eLeadClust.getEnergy() / primary.getEnergy() );
    helper.SetVariable( "fracSumBHCalVsBEMC", eSumECal / (eSumECal + eSumHCal) );
    helper.SetVariable( "fracLeadBHCalVsBEMC", eLeadClust.getEnergy() / (eLeadClust.getEnergy() + hLeadClust.getEnergy()) );
    helper.SetVariable( "diffSumBEMC", (eSumECal - primary.getEnergy()) / primary.getEnergy() );
    helper.SetVariable( "diffLeadBEMC", (eLeadClust.getEnergy() - primary.getEnergy()) / primary.getEnergy() );

    // if no energy in BHCal or BIC, skip event
    const bool isHCalNonzero = (eSumHCal > 0.);
    const bool isECalNonzero = (eSumECal > 0.);
    if (!isHCalNonzero && !isECalNonzero) continue;

    // ------------------------------------------------------------------------
    // scfi cluster/hit loops
    // ------------------------------------------------------------------------
    edm4eic::Cluster sLeadClust;

    // loop over scfi ecal clusters
    float eSumScFi  = 0.;
    float eLeadScFi = 0.;
    for (edm4eic::Cluster sClust : scfiClusters) {

      if (sClust.getEnergy() > eLeadScFi) {
        sLeadClust = sClust;
        eLeadScFi  = sClust.getEnergy();
      }
      eSumScFi += sClust.getEnergy();

    }  // end scfi cluster loop

    // fill scfi cluster variables
    helper.SetVariable( "nClustScFi", (float) scfiClusters.size() );
    helper.SetVariable( "eSumScFi", eSumScFi );
    helper.SetVariable( "eLeadScFi", sLeadClust.getEnergy() );
    helper.SetVariable( "hLeadScFi", edm4hep::utils::eta(sLeadClust.getPosition()) );
    helper.SetVariable( "fLeadScFi", edm4hep::utils::angleAzimuthal(sLeadClust.getPosition()) );

    // loop over scfi hits
    std::map<int32_t, float> mapScFiSumToLayer;
    for (edm4eic::CalorimeterHit sRecHit : scfiHits) {
      mapScFiSumToLayer[ sRecHit.getLayer() ] += sRecHit.getEnergy();
    }  // end scfi hit loop

    // fill scfi layer variables
    helper.SetVariable( "eSumScFiLayer1", mapScFiSumToLayer[1] );
    helper.SetVariable( "eSumScFiLayer2", mapScFiSumToLayer[2] );
    helper.SetVariable( "eSumScFiLayer3", mapScFiSumToLayer[3] );
    helper.SetVariable( "eSumScFiLayer4", mapScFiSumToLayer[4] );
    helper.SetVariable( "eSumScFiLayer5", mapScFiSumToLayer[5] );
    helper.SetVariable( "eSumScFiLayer6", mapScFiSumToLayer[6] );
    helper.SetVariable( "eSumScFiLayer7", mapScFiSumToLayer[7] );
    helper.SetVariable( "eSumScFiLayer8", mapScFiSumToLayer[8] );
    helper.SetVariable( "eSumScFiLayer9", mapScFiSumToLayer[9] );
    helper.SetVariable( "eSumScFiLayer10", mapScFiSumToLayer[10] );
    helper.SetVariable( "eSumScFiLayer11", mapScFiSumToLayer[11] );
    helper.SetVariable( "eSumScFiLayer12", mapScFiSumToLayer[12] );

    // ------------------------------------------------------------------------
    // imaging cluster loop
    // ------------------------------------------------------------------------
    edm4eic::Cluster iLeadClust;

    // loop over imaging ecal clusters (layers)
    float eSumImage  = 0.;
    float eLeadImage = 0.;
    for (edm4eic::Cluster iClust : imageClusters) {

      if (iClust.getEnergy() > eLeadImage) {
        iLeadClust = iClust;
        eLeadImage = iClust.getEnergy();
      }
      eSumImage += iClust.getEnergy();

    }  // end imaging cluster loop

    // fill imaging cluster variables
    helper.SetVariable( "nClustImage", (float) imageClusters.size() );
    helper.SetVariable( "eSumImage", eSumImage );
    helper.SetVariable( "eLeadImage", iLeadClust.getEnergy() );
    helper.SetVariable( "hLeadImage", edm4hep::utils::eta(iLeadClust.getPosition()) );
    helper.SetVariable( "fLeadImage", edm4hep::utils::angleAzimuthal(iLeadClust.getPosition()) );

    // loop over scfi hits
    std::map<int32_t, float> mapImageSumToLayer;
    for (edm4eic::CalorimeterHit iRecHit : imageHits) {
      mapImageSumToLayer[ iRecHit.getLayer() ] += iRecHit.getEnergy();
    }  // end scfi hit loop

    // fill image layer variables
    helper.SetVariable( "eSumImageLayer1", mapImageSumToLayer[1] );
    helper.SetVariable( "eSumImageLayer2", mapImageSumToLayer[2] );
    helper.SetVariable( "eSumImageLayer3", mapImageSumToLayer[3] );
    helper.SetVariable( "eSumImageLayer4", mapImageSumToLayer[4] );
    helper.SetVariable( "eSumImageLayer5", mapImageSumToLayer[5] );
    helper.SetVariable( "eSumImageLayer6", mapImageSumToLayer[6] );

    // ------------------------------------------------------------------------
    // fill ntuple
    //  -----------------------------------------------------------------------
    ntForCalib -> Fill( helper.GetValues().data() );

  }  // end frame loop
  std::cout << "    Finished frame loop" << std::endl;

  // save output & close files
  output     -> cd();
  ntForCalib -> Write();
  output     -> Close();

  // announce end & exit
  std::cout << "  End of macro!\n" << std::endl;
  return;

}

// end ========================================================================
