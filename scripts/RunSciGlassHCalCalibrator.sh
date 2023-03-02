#!/bin/bash
# -----------------------------------------------------------------------------
# 'RunSciGlassHCalCalibrator.sh'
# Derek Anderson
# 03.02.2023
#
# A simple script to run the
# JCalibrateHCal JANA plugin
# -----------------------------------------------------------------------------

# i/o files
input="../forECalStudy.sciGlass_run0.e5th35145n25Kpim.d22m2y2023.edm4hep.root"
podio="test_sciglass.podio.root"
output="test_sciglass.hists.root"

# output collections from EICrecon
collections="HcalBarrelRecHits,HcalBarrelClusters,EcalBarrelSciGlassClusters,HcalBarrelIslandProtoClusters,HcalBarrelTruthClusters,HcalBarrelTruthProtoClusters,GeneratedParticles"

# set geometry to SciGlass bemc
export DETECTOR_CONFIG=epic

eicrecon -Pplugins=JCalibrateHCalWithSciGlass -Ppodio:output_include_collections=$collections -Ppodio:output_file=$podio -Phistsfile=$output $input

# end -------------------------------------------------------------------------
