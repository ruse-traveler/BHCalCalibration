#!/bin/bash
# 'RunHCalCalibrator.sh'
# Derek Anderson
# 11.07.2022
#
# A simple script to run the
# JCalibrateHCal JANA plugin

# i/o files
input="../forVaryingClusterPars.e10th70n10Kpip.d19m1y2023.edm4hep.root"
podio="forPodioReaderTest_fromEicRecon.e5th70n500pip.d18m1y2023.podio.root"
output="forPodioReaderTest_fromEicRecon.e5th70n500pip.d18m1y2023.jcalibratehcal.root"

# output collections from EICrecon
collections="HcalBarrelRecHits,HcalBarrelClusters,HcalBarrelIslandProtoClusters,HcalBarrelTruthClusters,HcalBarrelTruthProtoClusters,GeneratedParticles"

# clustering parameters
scale="5,5"
distXY="150,150"
sector="5"
center="0.03"
energy="0.003"

eicrecon -Pplugins=JCalibrateHCal -PHCAL:HcalBarrelIslandProtoClusters:dimScaledLocalDistXY=$scale -PHCAL:HcalBarrelIslandProtoClsuters:localDistXY=$distXY -PHCAL:HcalBarrelIslandProtoClusters:sectorDist=$sector -PHCAL:HcalBarrelIslandProtoClusters:minClusterCenterEdep=$center -PHCAL:HcalBarrelIslandProtoClusters:minClusterHitEdep=$energy -Ppodio:output_include_collections=$collections -Ppodio:output_file=$podio -Phistsfile=$output $input

# end -------------------------------------------------------------------------
