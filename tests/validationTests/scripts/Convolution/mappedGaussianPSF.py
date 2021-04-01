import platosim.referenceFrames as rf
import numpy                    as np

from test                import Test
from math                import degrees, pow, sqrt
from platosim.validation import fitGaussian2D, gaussian2D
import matplotlib.pyplot as plt



"""This test checks checks three things:
1. It checks the simulation with a mapped gaussian PSF model. This is done by generating a theoretical sub pixel
   image and subtracting this from the sub pixel image obtained from the simulator.

2. It check that the charge diffusion spreads the image out. This is done by fitting the sub pixel image of a simulation with a large diffusion to a gaussian function. The test passes if the fitted standard deviation is much higher then the one we would expect from a simulation without charge diffusion.

3. It checks that the jitter smearing spreads the image out. This is done by comparing the pixel image with and without jitter smearing. 
"""





class MappedGaussianPSF(Test):

    def setNr(self):

        self.nr = "008.1"

    def setAllEffects(self):

        super().setAllEffects()
        self.sim["ObservingParameters/NumExposures"] = 1
        self.sim["CCD/IncludeConvolution"]           = "yes"

        self.sim["ControlHDF5Content/WriteSubPixelImages"] = "yes"
        self.sim["PSF/Model"] = "MappedGaussian"

        self.sim["SubField/NumRows"]    = 10
        self.sim["SubField/NumColumns"] = 10


        # Create a SkyMap with a star that would fall onto that point in the sub field.
        nRows     = self.sim["SubField/NumRows"]
        nColumns  = self.sim["SubField/NumColumns"]

        self.pRow = nRows / 2
        self.pCol = nColumns / 2

        ra, dec = rf.pixelToSkyCoordinates(self.sim, "2", self.pRow, self.pCol)

        ra  = degrees(ra)
        dec = degrees(dec)


        starCatalogFilename = self.outputDir + "/starCatalog"+ self.nr + ".txt"
        myFile = open(starCatalogFilename, "w")
        myFile.write("# RA DEC Vmag starID\n")
        myFile.write("{0}  {1}  {2}  {3}\n".format(ra, dec, 16.5, 1))
        myFile.close()

        self.sim["ObservingParameters/StarCatalogFile"] = starCatalogFilename



    def runSimulation(self):

        simFile = self.sim.run(removeOutputFile = True)
        self.image = simFile.getSubPixelImage(0)

        self.sim["PSF/MappedGaussian/IncludeChargeDiffusion"] = 'yes'
        self.sim["PSF/MappedGaussian/ChargeDiffusionStrength"] = 2
        simFile = self.sim.run(removeOutputFile = True)
        self.image2 = simFile.getSubPixelImage(0)

        self.sim["PSF/MappedGaussian/IncludeChargeDiffusion"] = 'no'
        self.sim["PSF/MappedGaussian/IncludeJitterSmoothing"] = "yes"
        simFile = self.sim.run(removeOutputFile = True)
        self.image3 = simFile.getSubPixelImage(0)


        





    def compareMapped(self):
        # We check that the image matches the predicted Gaussian distribution.

        subPixels      = self.sim["SubField/SubPixels"]
        colMax, rowMax = self.pCol * subPixels, self.pRow * subPixels
        amplitude      = np.max(self.image) - np.min(self.image)
        sigma          = self.sim["PSF/MappedGaussian/Sigma"] * subPixels
        function       = gaussian2D(amplitude, colMax, rowMax, sigma, sigma)
        columns, rows  = self.image.shape
        theoPrediction = np.array([[function(row, col) + np.min(self.image) for col in range(columns)] for row in range(rows)])

        difference     = theoPrediction - self.image


        return np.max(difference) / amplitude < 0.5



    
    def compareDiffusion(self):
        # We check that the image with diffusion strength 2, gives a Gaussian function with a
        # much higher standard deviation then the image without diffusion.

        subPixels = self.sim["SubField/SubPixels"]
        amplitude = np.max(self.image2) - np.min(self.image2)
        sigma0    = self.sim["PSF/MappedGaussian/Sigma"] * subPixels
        param = fitGaussian2D(self.image2, amplitude, self.pRow, self.pCol, sigma0, sigma0)
        return (param[-1] > 1.7 * sigma0) and (param[-2] > 1.7 * sigma0)



    
    def compareJitterSmooting(self):
        # We checks that the image with jitter smoothing has a higher standard deviation then
        # the image without jitter smoothing.
        
        subPixels  = self.sim["SubField/SubPixels"]
        sigmaPSFt  = self.sim["PSF/MappedGaussian/Sigma"]
               
        amplitude1 = np.max(self.image) -  np.min(self.image)
        amplitude3 = np.max(self.image3) - np.min(self.image3)

        sigma0     = sigmaPSFt * subPixels
        sigmaJ     = sqrt( pow(sigmaPSFt, 2) + pow( 0.5 / subPixels, 2)) * subPixels
        
        param  = fitGaussian2D(self.image,  amplitude1, self.pRow, self.pCol, sigma0, sigma0)
        param2 = fitGaussian2D(self.image3, amplitude3, self.pRow, self.pCol, sigmaJ, sigmaJ)

       
        param  = np.array(param[-1:-3:-1])
        param2 = np.array(param2[-1:-3:-1])

        difference = (param2 - param) / (sigmaJ - sigma0)
        condition  = np.logical_and(0.7 < difference, difference < 1.3)
        return np.all(condition)

       

        



    def compare(self):

        test1 = self.compareMapped()
        test2 = self.compareDiffusion()
        test3 = self.compareJitterSmooting()

        return test1 and test2 and test3






if __name__ == "__main__":
    t = MappedGaussianPSF()
    print(t.run())
