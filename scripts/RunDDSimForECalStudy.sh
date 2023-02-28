#!/bin/bash
# -----------------------------------------------------------------------------
# 'RunDDSimForECalStudy.sh'
# Derek Anderson
# 02.21.2023
#
# A simple script to run single particle
# events via dd4hep
# -----------------------------------------------------------------------------

# input parameters
sciglass='$DETECTOR_PATH/epic.xml'
imaging='$DETECTOR_PATH/epic_imaging.xml'
steerer="../steering.forECalStudy_e2th35145pim.py"

# output parameters
numEvt=25000
outIma="forECalStudy.imaging_run0.e2th35145n25Kpim.d21m2y2023.edm4hep.root"
outSci="forECalStudy.sciGlass_run0.e2th35145n25Kpim.d21m2y2023.edm4hep.root"

# run dd4hep
ddsim --steeringFile $steerer --compactFile $imaging -G -N $numEvt --outputFile $outIma
ddsim --steeringFile $steerer --compactFile $sciglass -G -N $numEvt --outputFile $outSci

# end -------------------------------------------------------------------------
