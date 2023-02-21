This repository contains code to evaluate the response of the simulated Hadronic Calorimeter of the ePIC calorimeter. There are a few pieces of code:

  - `JCalibrateHCal`: A simple JANA plugin to compare the reconstructed hit and cluster energy in the HCal to simulated particles.
  - `PCalibrateHCal.C`: A macro which uses PODIO's ROOTReader to perform the same analysis as `JCalibrateHCal` with the PODIO output from EICrecon

Both `JCalibrateHCal` and `PCalibrateHCal.C` must be run in the EIC environment.

### JCalibrateHCal Usage
```
# after compiling EICrecon, do:
eicmkplugin.py JCalibrateHCal
cp JCalibrateHCalProcessor.* $EICrecon_ROOT/JCalibrateHCal/
cmake -S JCalibrateHcal -B JCalibrateHCal/build
cmake --build JCalibrateHCal/build --target install
eicrecon -Pplugins=JCalibrateHCal <input edm4hep file>
```

### PCalibrateHCal Usage
```
root -b -q PCalibrateHCal.C
```

---

### TODO items:
  - [Minor] Refactor histogram creation/style setting
