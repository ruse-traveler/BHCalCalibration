from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV, degree
SIM = DD4hepSimulation()
SIM.gun.energy = 5*GeV
SIM.gun.particle = "pi+"
SIM.gun.position = (0.0, 0.0, 0.0)
SIM.gun.distribution = "cos(theta)"
SIM.gun.thetaMin = 70*degree
SIM.gun.thetaMax = 80*degree
