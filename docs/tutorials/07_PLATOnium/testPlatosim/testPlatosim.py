#!/usr/bin/env python3

import os
from platosim.simulation import Simulation

# Setup and run simulation
runDir  = os.getcwd()
sim     = Simulation('output', outputDir=runDir)
sim["ObservingParameters/NumExposures"] = 1
simFile = sim.run()
