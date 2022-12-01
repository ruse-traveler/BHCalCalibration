#!/bin/bash
# 'RunHCalCalibrator.sh'
# Derek Anderson
# 11.07.2022
#
# A simple script to run the
# JCalibrateHCal JANA plugin

input="../forHCalClusterCheck_p2d7080n5000pip_d12m11y2022.edm4hep.root"
collections="HcalBarrelRecHits,EcalBarrelSciGlassRecHits,HcalBarrelClusters,EcalBarrelSciGlassClusters,HcalBarrelTruthClusters,EcalBarrelSciGlassTruthClusters,GeneratedParticles"

eicrecon -Pplugins=JCalibrateHCal -Ppodio:output_include_collections=$collections $input
