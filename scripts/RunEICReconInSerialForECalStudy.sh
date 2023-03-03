#!/bin/bash
# -----------------------------------------------------------------------------
# 'RunEICReconInSerialForECalStudy.sh'
# Derek Anderson
# 03.03.20223
#
# A simple script to run EICrecon over
# a set of inputs. Requires 'RunSingle
# EICReconForECalStudy.sh'
#
# NOTE: this assumes EICrecon has been
# compiled beforehand.
# -----------------------------------------------------------------------------

# arrays for i/o files
declare -a inima
declare -a insci
declare -a podima
declare -a podsci
declare -a plugima
declare -a plugsci

# environment parameters
rundir="/sphenix/u/danderson/eic/EICrecon_test"
shellScript="/sphenix/u/danderson/scripts/eic-shell"
setupScript="/sphenix/u/danderson/eic/EICrecon_test/bin/eicrecon-this.sh"
installDir="/sphenix/u/danderson/eic/EICrecon_MY"

# detector parameters
detector="/sphenix/u/danderson/eic/epic/install/setup.sh"
sciglass="epic"
imaging="epic_imaging"
pluginSci="JCalibrateHCalWithSciGlass"
pluginIma="JCalibrateHCalWithImaging"
collectSci="HcalBarrelRecHits,HcalBarrelClusters,EcalBarrelSciGlassClusters,HcalBarrelIslandProtoClusters,HcalBarrelTruthClusters,HcalBarrelTruthProtoClusters,GeneratedParticles"
collectIma="HcalBarrelRecHits,HcalBarrelClusters,EcalBarrelImagingMergedClusters,HcalBarrelIslandProtoClusters,HcalBarrelTruthClusters,HcalBarrelTruthProtoClusters,GeneratedParticles"

# i/o parameters
indirIma="/sphenix/user/danderson/eic/ecal_study/ddsim_output/imaging"
indirSci="/sphenix/user/danderson/eic/ecal_study/ddsim_output/sciglass"
outdir="/sphenix/user/danderson/eic/ecal_study/eicrecon_output"

# input files
inima[0]="forECalStudy.imaging_run0.e5th35145n25Kpim.d22m2y2023.edm4hep.root"
inima[1]="forECalStudy.imaging_run1.e5th35145n25Kpim.d22m2y2023.edm4hep.root"
inima[2]="forECalStudy.imaging_run2.e5th35145n25Kpim.d22m2y2023.edm4hep.root"
inima[3]="forECalStudy.imaging_run3.e5th35145n25Kpim.d22m2y2023.edm4hep.root"

insci[0]="forECalStudy.sciGlass_run0.e5th35145n25Kpim.d22m2y2023.edm4hep.root"
insci[1]="forECalStudy.sciGlass_run1.e5th35145n25Kpim.d22m2y2023.edm4hep.root"
insci[2]="forECalStudy.sciGlass_run2.e5th35145n25Kpim.d22m2y2023.edm4hep.root"
insci[3]="forECalStudy.sciGlass_run3.e5th35145n25Kpim.d22m2y2023.edm4hep.root"

# output podio files
podima[0]="forECalStudy.imaging_run0.e5th35145n25Kpim.d2m3y2023.podio.root"
podima[1]="forECalStudy.imaging_run1.e5th35145n25Kpim.d2m3y2023.podio.root"
podima[2]="forECalStudy.imaging_run2.e5th35145n25Kpim.d2m3y2023.podio.root"
podima[3]="forECalStudy.imaging_run3.e5th35145n25Kpim.d2m3y2023.podio.root"

podsci[0]="forECalStudy.sciGlass_run0.e5th35145n25Kpim.d2m3y2023.podio.root"
podsci[1]="forECalStudy.sciGlass_run1.e5th35145n25Kpim.d2m3y2023.podio.root"
podsci[2]="forECalStudy.sciGlass_run2.e5th35145n25Kpim.d2m3y2023.podio.root"
podsci[3]="forECalStudy.sciGlass_run3.e5th35145n25Kpim.d2m3y2023.podio.root"

# output 
plugima[0]="forECalStudy.imaging_run0.e5th35145n25Kpim.d2m3y2023.plugin.root"
plugima[1]="forECalStudy.imaging_run1.e5th35145n25Kpim.d2m3y2023.plugin.root"
plugima[2]="forECalStudy.imaging_run2.e5th35145n25Kpim.d2m3y2023.plugin.root"
plugima[3]="forECalStudy.imaging_run3.e5th35145n25Kpim.d2m3y2023.plugin.root"

plugsci[0]="forECalStudy.sciGlass_run0.e5th35145n25Kpim.d2m3y2023.plugin.root"
plugsci[1]="forECalStudy.sciGlass_run1.e5th35145n25Kpim.d2m3y2023.plugin.root"
plugsci[2]="forECalStudy.sciGlass_run2.e5th35145n25Kpim.d2m3y2023.plugin.root"
plugsci[3]="forECalStudy.sciGlass_run3.e5th35145n25Kpim.d2m3y2023.plugin.root"

# loop over files
(( iFile=0 ))
for input in ${inima[@]}; do

  # create imaging driver script
  touch $rundir/DoImagingEICRecon.sh
  echo "#!/bin/bash" > $rundir/DoImagingEICRecon.sh
  echo "source $rundir/RunSingleEICReconForECalStudy.sh $setupScript $rundir $detector $imaging $installDir $collectIma $pluginIma $input ${podima[$iFile]} ${plugima[$iFile]}" > $rundir/DoImagingEICRecon.sh
  chmod u+x $rundir/DoImagingEICRecon.sh

  # copy imaging input to run directory
  cp $indirIma/$input $rundir

  # run eicrecon for imaging
  $shellScript -- $rundir/DoImagingEICRecon.sh
  mv $rundir/${podima[$iFile]} $outdir
  mv $rundir/${plugima[$iFile]} $outdir
  rm $input

  # create SciGlass driver script
  touch $rundir/DoSciGlassEICRecon.sh
  echo "#!/bin/bash" > $rundir/DoSciGlassEICRecon.sh
  echo "source $rundir/RunSingleEICReconForECalStudy.sh $setupScript $rundir $detector $sciglass $installDir $collectSci $pluginSci ${insci[$iFile]} ${podsci[$iFile]} ${plugsci[$iFile]}" > $rundir/DoSciGlassEICRecon.sh
  chmod u+x $rundir/DoSciGlassEICRecon.sh

  # copy SciGlass input to run directory
  cp $indirSci/${insci[$iFile]} $rundir

  # run eicrecon for SciGlass
  $shellScript -- $rundir/DoSciGlassEICRecon.sh
  mv $rundir/${podsci[$iFile]} $outdir
  mv $rundir/${plugsci[$iFile]} $outdir
  rm ${insci[$iFile]}

  # clean up & increment counter
  rm $rundir/DoImagingEICRecon.sh
  rm $rundir/DoSciGlassEICRecon.sh
  (( iFile++ ))

done

unset inima
unset insci
unset podima
unset podsci
unset plugima
unset plugsci

# end -------------------------------------------------------------------------
