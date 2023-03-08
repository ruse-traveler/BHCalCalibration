#!/bin/bash
# -----------------------------------------------------------------------------
# 'RunSingleEICReconForECalStudy.sh'
# Derek Anderson
# 03.03.2023
#
# Runs EICrecon for provided arguments.
# Called by 'RunEICReconInSerialForECal
# Study.sh'
#
# NOTE: this assumes EICrecon has been
# compiled beforehand.
# -----------------------------------------------------------------------------

# passed arguments
setup=$1
rundir=$2
detector=$3
config=$4
installdir=$5
collections=$6
plugin=$7
input=$8
podout=$9
plugout=${10}

# set eicrecon
source $setup
export PATH=$PATH':/usr/local/bin'

# set geometry
source $detector
export DETECTOR_CONFIG=$config

# set install path
export EICrecon_MY=$installdir
cd $rundir

# build plugin
cmake -S $plugin -B $plugin/build
cmake --build $plugin/build --target install

# run eicrecon
eicrecon -Pplugins=$plugin -Ppodio:output_include_collections=$collections -Ppodio:output_file=$podout -Phistsfile=$plugout $input
exit

# end -------------------------------------------------------------------------
