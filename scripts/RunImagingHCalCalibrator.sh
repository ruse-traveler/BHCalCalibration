#!/bin/bash
# -----------------------------------------------------------------------------
# 'RunImagingHCalCalibrator.sh'
# Derek Anderson
# 03.02.2023
#
# A simple script to run the
# JCalibrateHCal JANA plugin
# -----------------------------------------------------------------------------

# i/o files
input="../forECalStudy.imaging_run0.e5th35145n25Kpim.d22m2y2023.edm4hep.root"
podio="test_imaging.podio.root"
output="test_imaging.jcalibratehcal.root"

# output collections from EICrecon
collections="HcalBarrelRecHits,HcalBarrelClusters,EcalBarrelImagingMergedClusters,HcalBarrelIslandProtoClusters,HcalBarrelTruthClusters,HcalBarrelTruthProtoClusters,GeneratedParticles"

# set geometry to imaging bemc
export DETECTOR_CONFIG=epic_imaging

eicrecon -Pplugins=JCalibrateHCalWithImaging -Ppodio:output_include_collections=$collections -Ppodio:output_file=$podio -Phistsfile=$output $input

# end -------------------------------------------------------------------------
