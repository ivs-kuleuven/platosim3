#!/usr/bin/env python3

"""
This script illustrates a usage example for PlatoSim3 with the worker framework: here we want to use worker out-of-box in a so-called "parameter sweep", meaining that we want to run the same software (PlatoSim3) with a large amount of different input parameters (inputfils.yaml). As an example lets say we want to make substantial timeseries for each pointing of PLATO original lifetime. This case will be very efficient in the worker framework, if each timeseries takes a few minutes to an hour. If the number of pointings are sifficiently low, with 10 planed for PLATO initially, we can split each pointing-timeseries calculation to an individual core on the same node. To conclude, a text file with the 10 pointing coordinates (RA, Dec) will be used as input parameters to this python script, which will be executed by our job script.    
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

import setupEnviroment as env
from simulation import Simulation

# SETUP SIMULATION

# Input arguments
jobName     = sys.argv[0][:-3]
raPointing  = sys.argv[1]
decPointing = sys.argv[2]
inputFile   = env.projectDir + '/inputfiles/inputfile.yaml'
statCatalog = 1

# Setup simulation
sim = Simulation(jobName, configurationFile=inputFile, outputDir=env.workDir)

sim["ObservingParameters/RApointing"]   = raPointing
sim["ObservingParameters/DecPointing"]  = decPointing
sim["ObservingParameters/NumExposures"] = 1

# Run simulation
simFile = sim.run(removeOutputFile=True)
# img = simFile.getImage(0)
# plt.imshow(img, interpolation='nearest', origin='lower', cmap='cool')
# plt.show()
