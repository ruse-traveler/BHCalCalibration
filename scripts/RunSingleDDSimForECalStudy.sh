#!/bin/bash
# -----------------------------------------------------------------------------
# 'RunSingleDDSimForECalStudy.sh'
# Derek Anderson
# 02.24.2023
#
# Runs ddsim for provided arguments.
# Called by 'RunDDsimInSerialForEcalStudy.sh'
# -----------------------------------------------------------------------------

# passed arguments
setup=$1
rundir=$2
steerer=$3
compact=$4
numevts=$5
output=$6

# set up detector and go to run directory
source $setup
cd $rundir

# run dd4hep
ddsim --steeringFile $steerer --compactFile $compact -G -N $numevts --outputFile $output

# end -------------------------------------------------------------------------
