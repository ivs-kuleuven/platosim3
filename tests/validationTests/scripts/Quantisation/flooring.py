from test import Test
import numpy as np



class Flooring(Test):
    """This test is designed to check the flooring effect in PlatoSim. The test runs the simulation with and without quantisation. 
    The two obtained images are then subtracted and all the values of the corresponding image should be in the interval [0, -1[. 
    The values in the image with flooring should be of the type int and the image without flooring should be of the type float."""

    

    def setNr(self):
        self.nr = "020.3"

    def setAllEffects(self):

        super().setAllEffects()

        self.sim["ObservingParameters/NumExposures"] = 1
        self.sim["SubField/NumRows"]    = 4510
        self.sim["SubField/NumColumns"] = 4510

        self.sim["ObservingParameters/StarCatalogFile"] = self.inputDir + "/starcatalog.txt"
        

        # Reset the Camera so that the stars in the input file fall onto the CCD.
        self.sim["CCD/Position"] = "Custom"
        self.sim["Telescope/GroupID"] = "Custom"
        self.sim["Telescope/AzimuthAngle"] = 0.
        self.sim["Telescope/TiltAngle"]    = 0.
        
        # Enable quantisation (only flooring)
        self.sim["CCD/IncludeQuantisation"] = "yes"
        self.sim["FEE/Gain/RefValueLeft"]   = 1
        self.sim["FEE/Gain/RefValueRight"]  = 1
        self.sim['CCD/Gain/RefValueLeft']   = 1
        self.sim['CCD/Gain/RefValueRight']  = 1
        self.sim['FEE/ElectronicOffset/RefValue'] = 0

        self.sim['CCD/IncludeDigitalSaturation'] = 'yes'
        self.saturationLimit = self.sim['CCD/DigitalSaturation']

    def runSimulation(self):

        # Run the simulation with quantisation
        simFileWithFlooring = self.sim.run(removeOutputFile = True)
        self.imageWithFlooring = simFileWithFlooring.getImage(0)

        # Run the simulation without quantisation
        self.sim['CCD/IncludeQuantisation'] = 'no'
        simFileWithoutFlooring = self.sim.run(removeOutputFile = True)
        self.imageWithoutFlooring = simFileWithoutFlooring.getImage(0)

    def compare(self):

        imageFlooring   = self.imageWithFlooring
        imageNoFlooring = self.imageWithoutFlooring

        # We look at the difference between the two images at the pixels where the saturation limit is not reached
        difference = imageNoFlooring[imageNoFlooring <= self.saturationLimit] -  imageFlooring[imageNoFlooring <= self.saturationLimit]

        condition1 = np.max(difference) <= 1 and np.min(difference) >= 0
        condition2 = imageFlooring.dtype == np.uint16 and imageNoFlooring.dtype == np.float32
        
        return condition1 and condition2







if __name__ == "__main__":
    t = Flooring()
    print(t.run())
