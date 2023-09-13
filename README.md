This repository contains code to evaluate the response of the simulated ePIC Barrel Hadronic Calorimeter (BHCal). There are a few pieces of code:

  - `JCalibrateHCalWithImaging`: A simple JANA plugin to compare the reconstructed hit and cluster energy in the BHCal and BECal to simulated particles. This prepares a TNtuple to be read in by a ROOT macro to train a TMVA model.
  - `TrainAndApplyBHCalCalibration.cxx`: A ROOT macro which reads in the TNtuple from `JCalibrateHCalWithImaging` and trains a TMVA model according to set specifications.

### JCalibrateHCal Usage
```
# after compiling EICrecon, do:
eicmkplugin.py JCalibrateHCal
cp JCalibrateHCalProcessor.* $EICrecon_ROOT/JCalibrateHCal/
cmake -S JCalibrateHcal -B JCalibrateHCal/build
cmake --build JCalibrateHCal/build --target install
eicrecon -Pplugins=JCalibrateHCal <input edm4hep file>
```

### TrainAndApplyBHCalCalibration Usage
```
root -b -q PCalibrateHCal.C
```

---

### TODO items:
  - [Major] Streamline training/application workflow
  - [Major] Factor out resolution calculation in training/application macro
  - [Major] Factor out histogram plotting operations in trainining/application macro
  - [Minor] Refactor histogram creation/style setting
