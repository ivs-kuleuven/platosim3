from test import Test

import numpy as np






"""
This test is designed to check the field distortion. The test runs the simulation with and without polarization. 
We devide the obtained images and find their max and minimum values. These are compared the expected polarziation 
from the input file. The test passes: 
1. if the standard deviation is smaller then 0.01
2. The mean value is the one expected from the input file +- standard deviation.
"""




class Polarization(Test):

    def setNr(self):
        self.nr = "011.3"



        
    def setAllEffects(self):

        super().setAllEffects()
        self.sim["ObservingParameters/NumExposures"] = 1
        self.sim["SubField/NumRows"] = 2000
        self.sim["SubField/NumColumns"] = 2000



        
    def runSimulation(self):
        self.sim["CCD/IncludePolarization"] = "yes"
        self.simFile1 = self.sim.run(removeOutputFile = True)

        self.sim["CCD/IncludePolarization"] = "no"
        self.simFile2 = self.sim.run(removeOutputFile = True)




        
        
    def compare(self):
        quotient = self.simFile1.getImage(0) / self.simFile2.getImage(0)
        polarization = self.sim["CCD/Polarization/ExpectedValue"]
        
        minRatio = np.min(quotient)
        maxRatio = np.max(quotient)
        avgRatio = np.mean(quotient)
        stdRatio = np.std(quotient)

        condition1 = stdRatio < 0.01
        condition2 = minRatio - stdRatio < polarization < maxRatio + stdRatio
        
        
        print(condition1)
        print(condition2)
        return condition1 and condition2





    

if __name__ == "__main__":
    t = Polarization()
    print(t.run())
