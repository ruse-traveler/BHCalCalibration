#!/bin/bash
# 'RunHCalCalibrator.sh'
# Derek Anderson
# 11.07.2022
#
# A simple script to run the
# JCalibrateHCal JANA plugin

input="../hcal/input/nov_sim_campaign/pi-_2GeV_45to135deg.0001.eicrecon.tree.edm4eic.root"

eicrecon -Pplugins=JCalibrateHCal $input
