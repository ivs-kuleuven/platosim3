from test import Test
from platosim.validation import fitGaussian2D

import numpy as np

import matplotlib.pyplot as plt



"""
This test is designed to check the Brighter-Fatter Effect. The simulation is run for a star in the middle of the subfield with and without the BFE effect, for different 
magnitudes of the star. The test passes if the difference in width between the two simulations decreases with increasing magnitude, and the width is independent of the magnitude
for the simulations without the BFE effect.
"""











class BrighterFatterEffect(Test):

    def setNr(self):
        self.nr = "013"
    
    def setAllEffects(self):

        super().setAllEffects()
        self.sim["ObservingParameters/NumExposures"] = 1
        self.sim["CCD/IncludeConvolution"]           = "yes"

        self.dim                        = 9
        self.sim["SubField/NumRows"]    = self.dim
        self.sim["SubField/NumColumns"] = self.dim

        self.sim["ControlHDF5Content/WriteSubPixelImages"] = "yes"
        self.numSubPixels = self.sim["SubField/SubPixels"]

        self.sim["PSF/Model"] = "MappedGaussian"

    def setEffectsForNewTest(self):

        super().setAllEffects()
        self.sim["ObservingParameters/NumExposures"] = 1
        self.sim["CCD/IncludeConvolution"]           = "yes"

        self.sim["SubField/NumRows"]    = 1000
        self.sim["SubField/NumColumns"] = 1000

        self.sim["CCD/IncludePhotonNoise"] = "yes"
        self.sim["CCD/IncludeBFE"]         = "yes"

        # Don't want a star in the subfield.
        self.sim["ObservingParameters/DecPointing"] = - self.sim["ObservingParameters/DecPointing"]

    def runNewTest(self):

        backgrounds = np.arange(50, 2500, 100)
        means = np.array([])
        stds  = np.array([])

        for background in backgrounds:

            self.sim["Sky/SkyBackground"] = background
            simFile = self.sim.run(removeOutputFile = True)
            image = simFile.getImage(0)
            stds  = np.append(stds, np.std(image)**2)
            means = np.append(means, np.mean(image))
            #print(means)

        plt.plot(means, stds)
        plt.show()
        
            
        

        




    def runSimulation(self):
        
        magnitudes = np.arange(8, 20, 0.5)
        position   = np.array([self.dim / 2])
        
        starCatalogFilename = self.outputDir + "/starCatalog" + self.nr + ".txt"
        self.sigma               = self.sim["PSF/MappedGaussian/Sigma"]

        widthWithBFE     = np.array([]), np.array([])
        widthWithoutBFE  = np.array([]), np.array([])
        
        for mag in magnitudes:

            self.sim.createStarCatalogFileFromPixelCoordinates(position, position, np.array([mag]), np.array([1]), starCatalogFilename)

            # Run the simulation with BFE.
            self.sim["CCD/IncludeBFE"] = "yes"
            widthWithBFE = self.storeWidthAndTotalFlux(widthWithBFE)
            
            # Run the similation without BFE.
            self.sim["CCD/IncludeBFE"] = "no"
            widthWithoutBFE = self.storeWidthAndTotalFlux(widthWithoutBFE)

        self.widthWithBFE    = widthWithBFE
        self.widthWithoutBFE = widthWithoutBFE

        



    def storeWidthAndTotalFlux(self, width):

        sigma    = self.sigma
        simFile  = self.sim.run(removeOutputFile = True)
        image    = simFile.getImage(0)

        params   = fitGaussian2D(image, np.max(image), self.dim / 2, self.dim / 2, sigma, sigma, subtractConstant = True)
        width    = np.append(width, params[3])


        return width
        
        
            
    def compare(self):

        difference = self.widthWithBFE - self.widthWithoutBFE
        isFatter   = all(difference[:-1] > difference[1:])

        varWithout  = np.std(self.widthWithoutBFE)
        withoutIsConst = varWithout < 0.3

        
        return isFatter and withoutIsConst

        

       




if __name__ == "__main__":
    t = BrighterFatterEffect()
    #print(t.run())
    t.setEffectsForNewTest()
    t.runNewTest()
