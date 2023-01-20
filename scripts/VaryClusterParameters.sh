#!/bin/bash
# 'VaryClusterParameter.sh'
# Derek Anderson
# 01.19.2023
#
# This script runs the JCalibrateHCal
# JANA plugin while varying the
# EICrecon clustering parameters.

# declare lists
declare -a scale
declare -a center
declare -a energy
declare -a podio
declare -a output

# maximum dimension-scaled distance
scale[0]="5,5"
scale[1]="7,7"
scale[2]="10,10"
scale[3]="5,5"
scale[4]="7,7"
scale[5]="10,10"
scale[6]="5,5"
scale[7]="7,7"
scale[8]="10,10"
scale[9]="5,5"
scale[10]="7,7"
scale[11]="10,10"
scale[12]="5,5"
scale[13]="7,7"
scale[14]="10,10"
scale[15]="5,5"
scale[16]="7,7"
scale[17]="10,10"

# minimum central energy deposit
center[0]="0.03"
center[1]="0.03"
center[2]="0.03"
center[3]="0.05"
center[4]="0.05"
center[5]="0.05"
center[6]="0.1"
center[7]="0.1"
center[8]="0.1"
center[9]="0.03"
center[10]="0.03"
center[11]="0.03"
center[12]="0.05"
center[13]="0.05"
center[14]="0.05"
center[15]="0.1"
center[16]="0.1"
center[17]="0.1"

# minimum energy deposit
energy[0]="0.003"
energy[1]="0.003"
energy[2]="0.003"
energy[3]="0.003"
energy[4]="0.003"
energy[5]="0.003"
energy[6]="0.003"
energy[7]="0.003"
energy[8]="0.003"
energy[9]="0.05"
energy[10]="0.05"
energy[11]="0.05"
energy[12]="0.05"
energy[13]="0.05"
energy[14]="0.05"
energy[15]="0.05"
energy[16]="0.05"
energy[17]="0.05"

# podio output files
podio[0]="varyingClusterPars.dim5cen30adj3.e10th70n10Kpip.d19m1y2023.podio.root"
podio[1]="varyingClusterPars.dim7cen30adj3.e10th70n10Kpip.d19m1y2023.podio.root"
podio[2]="varyingClusterPars.dim10cen30adj3.e10th70n10Kpip.d19m1y2023.podio.root"
podio[3]="varyingClusterPars.dim5cen50adj3.e10th70n10Kpip.d19m1y2023.podio.root"
podio[4]="varyingClusterPars.dim7cen50adj3.e10th70n10Kpip.d19m1y2023.podio.root"
podio[5]="varyingClusterPars.dim10cen50adj3.e10th70n10Kpip.d19m1y2023.podio.root"
podio[6]="varyingClusterPars.dim5cen100adj3.e10th70n10Kpip.d19m1y2023.podio.root"
podio[7]="varyingClusterPars.dim7cen100adj3.e10th70n10Kpip.d19m1y2023.podio.root"
podio[8]="varyingClusterPars.dim10cen100adj3.e10th70n10Kpip.d19m1y2023.podio.root"
podio[9]="varyingClusterPars.dim5cen30adj50.e10th70n10Kpip.d19m1y2023.podio.root"
podio[10]="varyingClusterPars.dim7cen30adj50.e10th70n10Kpip.d19m1y2023.podio.root"
podio[11]="varyingClusterPars.dim10cen30adj50.e10th70n10Kpip.d19m1y2023.podio.root"
podio[12]="varyingClusterPars.dim5cen50adj50.e10th70n10Kpip.d19m1y2023.podio.root"
podio[13]="varyingClusterPars.dim7cen50adj50.e10th70n10Kpip.d19m1y2023.podio.root"
podio[14]="varyingClusterPars.dim10cen50adj50.e10th70n10Kpip.d19m1y2023.podio.root"
podio[15]="varyingClusterPars.dim5cen100adj50.e10th70n10Kpip.d19m1y2023.podio.root"
podio[16]="varyingClusterPars.dim7cen100adj50.e10th70n10Kpip.d19m1y2023.podio.root"
podio[17]="varyingClusterPars.dim10cen100adj50.e10th70n10Kpip.d19m1y2023.podio.root"

# JCalibrate HCal output files
output[0]="varyingClusterPars.dim5cen30adj3.e10th70n10Kpip.d19m1y2023.hists.root"
output[1]="varyingClusterPars.dim7cen30adj3.e10th70n10Kpip.d19m1y2023.hists.root"
output[2]="varyingClusterPars.dim10cen30adj3.e10th70n10Kpip.d19m1y2023.hists.root"
output[3]="varyingClusterPars.dim5cen50adj3.e10th70n10Kpip.d19m1y2023.hists.root"
output[4]="varyingClusterPars.dim7cen50adj3.e10th70n10Kpip.d19m1y2023.hists.root"
output[5]="varyingClusterPars.dim10cen50adj3.e10th70n10Kpip.d19m1y2023.hists.root"
output[6]="varyingClusterPars.dim5cen100adj3.e10th70n10Kpip.d19m1y2023.hists.root"
output[7]="varyingClusterPars.dim7cen100adj3.e10th70n10Kpip.d19m1y2023.hists.root"
output[8]="varyingClusterPars.dim10cen100adj3.e10th70n10Kpip.d19m1y2023.hists.root"
output[9]="varyingClusterPars.dim5cen30adj50.e10th70n10Kpip.d19m1y2023.hists.root"
output[10]="varyingClusterPars.dim7cen30adj50.e10th70n10Kpip.d19m1y2023.hists.root"
output[11]="varyingClusterPars.dim10cen30adj50.e10th70n10Kpip.d19m1y2023.hists.root"
output[12]="varyingClusterPars.dim5cen50adj50.e10th70n10Kpip.d19m1y2023.hists.root"
output[13]="varyingClusterPars.dim7cen50adj50.e10th70n10Kpip.d19m1y2023.hists.root"
output[14]="varyingClusterPars.dim10cen50adj50.e10th70n10Kpip.d19m1y2023.hists.root"
output[15]="varyingClusterPars.dim5cen100adj50.e10th70n10Kpip.d19m1y2023.hists.root"
output[16]="varyingClusterPars.dim7cen100adj50.e10th70n10Kpip.d19m1y2023.hists.root"
output[17]="varyingClusterPars.dim10cen100adj50.e10th70n10Kpip.d19m1y2023.hists.root"

# input file
input="../forVaryingClusterPars.e10th70n10Kpip.d19m1y2023.edm4hep.root"

# output collections from EICrecon
collections="HcalBarrelRecHits,HcalBarrelClusters,HcalBarrelIslandProtoClusters,HcalBarrelTruthClusters,HcalBarrelTruthProtoClusters,GeneratedParticles"

# loop over combinations
(( nCombo=0 ))
for combo in ${scale[@]}; do
  eicrecon -Pplugins=JCalibrateHCal -PHCAL:HcalBarrelIslandProtoClusters:dimScaledLocalDistXY=$combo -PHCAL:HcalBarrelIslandProtoClusters:minClusterCenterEdep=${center[$nCombo]} -PHCAL:HcalBarrelIslandProtoClusters:minClusterHitEdep=${energy[$nCombo]} -Ppodio:output_include_collections=$collections -Ppodio:output_file=${podio[$nCombo]} -Phistsfile=${output[$nCombo]} $input
  (( nCombo++ ))
done

# delete arrays
unset scale
unset center
unset energy
unset podio
unset output

# end -------------------------------------------------------------------------
