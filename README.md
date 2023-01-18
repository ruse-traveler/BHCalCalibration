This repository contains code to evaluate the response of the simulated Hadronic Calorimeter of the ePIC calorimeter. There are a few pieces of code:

  - `JCalibrateHCal`: A simple JANA plugin to compare the reconstructed hit and cluster energy in the HCal to simulated particles.

## JCalibrateHCal Usage
```
# after compiling EICrecon, do:
eicmkplugin.py JCalibrateHCal
cp JCalibrateHCalProcessor.* $EICrecon_ROOT/JCalibrateHCal/
cmake -S JCalibrateHcal -B JCalibrateHCal/build
cmake --build JCalibrateHCal/build --target install
eicrecon -Pplugins=JCalibrateHCal <input edm4hep file>
```

---

### TODO items:
  - [Major] Import functionality into ROOT macro to work with simulation campaing output
  - [Minor] Refactor histogram creation/style setting
