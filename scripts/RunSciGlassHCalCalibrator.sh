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
input="../forVaryingClusterPars.e10th70n10Kpip.d19m1y2023.edm4hep.root"
podio="forPodioReaderTest_fromEicRecon.e5th70n500pip.d18m1y2023.podio.root"
output="forPodioReaderTest_fromEicRecon.e5th70n500pip.d18m1y2023.jcalibratehcal.root"

# output collections from EICrecon
collections="HcalBarrelRecHits,HcalBarrelClusters,EcalBarrelSciGlassClusters,HcalBarrelIslandProtoClusters,HcalBarrelTruthClusters,HcalBarrelTruthProtoClusters,GeneratedParticles"

eicrecon -Pplugins=JCalibrateHCalWithSciGlass -Ppodio:output_include_collections=$collections -Ppodio:output_file=$podio -Phistsfile=$output $input

# end -------------------------------------------------------------------------
