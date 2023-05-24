#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
#from test import eprint
#from test       import Test

import platosim.referenceFrames as rf
from platosim.simulation import Simulation
#from platosim.validation import equatorial2galactic, galactic2equatorial, aberration


"""
This test checks the absolute aberration. The test consists out of 3 parts:

1. The first test plots the position of an aberrated star over one
year and a theoretical predicted path (where we assume a circular orbit with a constant velocity for the spacecraft) over one year. This figure is
saved into the ioFiles directory, it is possible to visually confirm that the two paths are similar.

2. The second test checks that the aberration is zero for a star whose position wrt the spacecraft is orthogonal to the spacecrafts velocity.

3. The last test checks that the aberration is maximal for a star whose position is parallel to the velocity.
"""


class DifferentialAberration():

    def simulate(self, DKA, pos, beginExposureNr, outputName):

        # Set up a simulation object
        sim = Simulation(outputName, outputDir=self.outputDir)

        # Define and set star catalog file
        starFileName = f'{self.outputDir}/starCatalog{self.nr}_{DKA}_pos{pos}_exp{beginExposureNr}.txt'
        row = col = np.array([pos])
        mag = np.array([11])
        sid = np.array([1])
        sim.createStarCatalogFileFromPixelCoordinates(row+4, col+4, mag, sid, starFileName)

        # Simulate only a small subfield
        sim["SubField/NumRows"]         = 8
        sim["SubField/NumColumns"]      = 8
        sim["SubField/ZeroPointRow"]    = pos
        sim["SubField/ZeroPointColumn"] = pos

        # Simulate a single exposure
        sim["ObservingParameters/NumExposures"]    = 1
        sim["ObservingParameters/BeginExposureNr"] = beginExposureNr
        
        # Activate kinematic differential aberration 
        sim["Camera/IncludeAberrationCorrection"] = DKA
        sim["Camera/AberrationCorrection/Type"]   = "differential"

        # Turn off all effect that can cause pixel displacements
        sim["Platform/UseJitter"] = "no"
        sim["Telescope/UseDrift"] = "no"

        # Run simulation
        f = sim.run(removeOutputFile=True)

        # Fetch parameters
        row, col = f.getStarPositions(sid)

        return np.sqrt((4-row)**2 + (4-col)**2)

        
        
    def run(self):

        self.nr = "003.2"
        self.outputDir = os.environ["PLATO_PROJECT_HOME"] + "/tests/validationTests/ioFiles/test" + self.nr
        tsep = int(50*86400/25)
        tend = int(360*86400/25)
        exp = np.arange(0,tend,tsep)
        pos = np.arange(10,3500,500)
        DKA = "yes"

        npos = len(pos)
        nexp = len(exp)
        c = np.zeros((npos,nexp))

        for i in range(npos):
            for j in range(nexp):
                print(i)
                c[i,j] = self.simulate(DKA, pos[i], exp[j], f"output{self.nr}")

        # Plot results
        plt.figure(figsize=(10,8))
        for i in range(npos):
            plt.plot(exp*25/86400, c[i,:], 'o-', label=f'(row,col) = ({pos[i]}, {pos[i]}) pix')
        plt.xlabel('Time [days]')
        plt.ylabel('Delta position [pixel]')
        plt.legend()
        plt.show()


if __name__ == "__main__":
    t = DifferentialAberration()
    print(t.run())

