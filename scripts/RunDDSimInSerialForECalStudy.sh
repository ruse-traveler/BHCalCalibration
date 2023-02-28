#!/bin/bash
# -----------------------------------------------------------------------------
# 'RunDDSimInSerialForECalStudy.sh'
# Derek Anderson
# 02.24.2023
#
# A simple script to run single particle
# events via dd4hep. Requires 'RunSingle
# DDsimForECalStudy.sh'
# -----------------------------------------------------------------------------

# arrays for output files
declare -a outima
declare -a outsci

# environment parameters
rundir="/sphenix/u/danderson/eic/run_ele"
shellScript="/sphenix/u/danderson/scripts/eic-shell"
setupScript="/sphenix/u/danderson/eic/epic/install/setup.sh"

# input parameters
sciglass="/sphenix/u/danderson/eic/epic/install/share/epic/epic.xml"
imaging="/sphenix/u/danderson/eic/epic/install/share/epic/epic_imaging.xml"
steerer="/sphenix/u/danderson/eic/steering.forECalStudy_e2th35145ele.py"

# output parameters
numevt=10000
outdir="/sphenix/user/danderson/eic/ecal_study/ddsim_output"

outima[0]="forECalStudy.imaging_run0.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outima[1]="forECalStudy.imaging_run1.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outima[2]="forECalStudy.imaging_run2.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outima[3]="forECalStudy.imaging_run3.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outima[4]="forECalStudy.imaging_run4.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outima[5]="forECalStudy.imaging_run5.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outima[6]="forECalStudy.imaging_run6.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outima[7]="forECalStudy.imaging_run7.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outima[8]="forECalStudy.imaging_run8.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outima[9]="forECalStudy.imaging_run9.e2th35145n10Kele.d27m2y2023.edm4hep.root"

outsci[0]="forECalStudy.sciGlass_run0.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outsci[1]="forECalStudy.sciGlass_run1.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outsci[2]="forECalStudy.sciGlass_run2.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outsci[3]="forECalStudy.sciGlass_run3.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outsci[4]="forECalStudy.sciGlass_run4.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outsci[5]="forECalStudy.sciGlass_run5.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outsci[6]="forECalStudy.sciGlass_run6.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outsci[7]="forECalStudy.sciGlass_run7.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outsci[8]="forECalStudy.sciGlass_run8.e2th35145n10Kele.d27m2y2023.edm4hep.root"
outsci[9]="forECalStudy.sciGlass_run9.e2th35145n10Kele.d27m2y2023.edm4hep.root"


# loop over files
(( iRun=0 ))
for output in ${outima[@]}; do

  # create imaging driver script
  touch $rundir/DoImagingDDSim.sh
  echo "#!/bin/bash" > $rundir/DoImagingDDSim.sh
  echo "source $rundir/RunSingleDDSimForECalStudy.sh $setupScript $rundir $steerer $imaging $numevt $output" > $rundir/DoImagingDDSim.sh
  chmod u+x $rundir/DoImagingDDSim.sh

  # run ddsim for imaging
  $shellScript -- $rundir/DoImagingDDSim.sh
  mv $rundir/$output $outdir

  # create imaging driver script
  touch $rundir/DoSciGlassDDSim.sh
  echo "#!/bin/bash" > $rundir/DoSciGlassDDSim.sh
  echo "source $rundir/RunSingleDDSimForECalStudy.sh $setupScript $rundir $steerer $sciglass $numevt ${outsci[$iRun]}" > $rundir/DoSciGlassDDSim.sh
  chmod u+x $rundir/DoSciGlassDDSim.sh

  # run dd4hep for SciGlass
  $shellScript -- $rundir/DoSciGlassDDSim.sh
  mv $rundir/${outsci[$iRun]} $outdir

  # clean up & increment counter
  rm $rundir/DoImagingDDSim.sh
  rm $rundir/DoSciGlassDDSim.sh
  (( iRun++ ))

done

# delete arrays
unset outima
unset outsci

# end -------------------------------------------------------------------------
