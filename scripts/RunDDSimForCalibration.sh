#!/bin/bash
# 'RunDDSimForCalibration.sh'
# Derek Anderson
# 12.01.2022
#
# A simple script to run dd4hep to generate
# an edm4hep file for JCalibrateHCal to
# work with

# input parameters
steerer="../steering.forHCalClusterCheck_p2th70pip.py"
compact='$DETECTOR_PATH/$DETECTOR_CONFIG.xml'

# output parameters
numEvt=10
output="test.edm4hep.root"

# run dd4hep
ddsim --steeringFile $steerer --compactFile $compact -G -N $numEvt --outputFile $output

# end -------------------------------------------------------------------------
