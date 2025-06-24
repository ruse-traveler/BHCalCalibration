**Note:** this is a public archive. This code is now maintained as part of [ruse-traveler/EpicBHCalPTDRStudies](https://github.com/ruse-traveler/EpicBHCalPTDRStudies) and [eic/snippets](https://github.com/eic/snippets).

This repository contains code to evaluate the response of the simulated ePIC Barrel Hadronic Calorimeter (BHCal).  There are a few pieces of code:

  - `JCalibrateHCalWithImaging`: Found in the `plugin` directory, this is a simple JANA plugin to compare the reconstructed hit and cluster energy in the BHCal and BECal to simulated particles. This also prepares a TNtuple to be read in by a ROOT macro to train a TMVA model.
  - `BHCalCalibration:` The source code is found in the `src` directory which is run with `DoBHCalCalibration.cxx`.  This ingests the TNtuple from `JCalibrateHCalWithImaging` and trains a TMVA model according to set specifications.

There is also `TrainAndApplyBHCalCalibration.cxx`, a ROOT macro which does the same thing as the `BHCalCalibration` package.  It was put together to test training and applying a TMVA model in the same macro.

The older macros which trained and applied the TMVA model separately (`DoHCalCalibration.C` and `ApplyHCalCalibration.C` respectively) are currently kept the `macros` directory for posterity.

### JCalibrateHCal Usage

After compiling `EICrecon`, create and compile the plugin with:

```
eicmkplugin.py JCalibrateHCal
cp <path to this repo>/plugin/JCalibrateHCalProcessor.* ./JCalibrateHCal/
cmake -S JCalibrateHcal -B JCalibrateHCal/build
cmake --build JCalibrateHCal/build --target install
```

Then run it with:

```
eicrecon -Pplugins=JCalibrateHCal <input edm4hep file>
```

### BHCalCalibration Usage

Compile the source code first with:

```
cd src/
./root-build
```

And then run the code with:

```
root -b -q DoBHCalCalibration.cxx
```

### TrainAndApplyBHCalCalibration Usage

No need to compile this beforehand, just run it in the usual manner.

```
root -b -q TrainAndApplyBHCalCalibration.cxx
```
