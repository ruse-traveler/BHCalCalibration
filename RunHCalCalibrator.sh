#!/bin/bash
# 'RunHCalCalibrator.sh'
# Derek Anderson
# 11.07.2022
#
# A simple script to run the
# JCalibrateHCal JANA plugin

input="../hcal/input/edm4hep/ohCalTest.p2t5n5000pip.d31m10y2022.edm4hep.root"
collections="HcalBarrelRawHits,HcalBarrelRecHits,HcalBarrelClusters,GeneratedParticles"

eicrecon -Pplugins=JCalibrateHCal -Ppodio:output_include_collections=$collections $input
