import numpy as np
from test       import Test
from platosim.validation import fitGaussian2D







"""
This test is designed to test 3 aspects of the simulation with an analytical non-Gaussian PSF. It tests that:
1. the spread of the star increases with respect to the distance to the optical axis, when sigma is constant. 
2. the behavour of the spread of the start with respect to time when sigma is given by an inputfile. 
3. it compares two simulations with and without charge diffusion.
"""









class AnalyticNonGaussianPSF(Test):



    
    def setNr(self):

        self.nr = "008.2"




        
    def setAllEffects(self):
        # Set the effects for the first simulation. 
        
        super().setAllEffects()
        self.sim["ObservingParameters/NumExposures"]     = 1
        self.sim["CCD/IncludeConvolution"]               = "yes"
        self.sim["PSF/Model"]                            = "AnalyticNonGaussian"
        self.sim["PSF/AnalyticNonGaussian/Sigma/Source"] = "ConstantValue"


        self.sim["SubField/NumRows"]    = 4510
        self.sim["SubField/NumColumns"] = 4510

        self.sim["PSF/AnalyticNonGaussian/Sigma/ConstantValue"] = 50



        # Create a star catalog file with multiple stars aligned over the diagonal of the subfield. 
        starCatalogFilename = self.outputDir + "/starCatalogFile1"+ self.nr + ".txt"
        
        self.xPosition1 = np.array([ 600, 1305, 2010])
        self.yPosition1 = np.array([3910, 3205, 2500])
                
        xPosition = self.xPosition1
        yPosition = self.yPosition1
        
        magnitude = 12.5 * np.ones(len(xPosition))
        index     = np.arange(1, len(xPosition) + 1)

        self.sim.createStarCatalogFileFromPixelCoordinates(yPosition, xPosition, magnitude, index, starCatalogFilename)
        self.sim["ObservingParameters/StarCatalogFile"] = starCatalogFilename
        




        



    def setEffectForPSFFromFile(self):
        # Sets the effects for the second simulation. A PSF inputfile is generated as well.
        
        self.sim["ObservingParameters/NumExposures"]     = 20
        self.sim["PSF/AnalyticNonGaussian/Sigma/Source"] = "FromFile"
        self.sim["SubField/NumRows"]    = 120
        self.sim["SubField/NumColumns"] = 120


        # Make star catalog for one star in the middle of the subfield.
        starCatalogFilename = self.outputDir + "/starCatalogFile2"+ self.nr + ".txt"
        
        self.xPosition2 = np.array([60])
        self.yPosition2 = np.array([60])
                
        xPosition = self.xPosition2
        yPosition = self.yPosition2
        
        magnitude = np.array([12.5])
        index     = np.array([1])

        self.sim.createStarCatalogFileFromPixelCoordinates(yPosition, xPosition, magnitude, index, starCatalogFilename)
        self.sim["ObservingParameters/StarCatalogFile"] = starCatalogFilename

        # Create an input file for the sigmas.
        psfInputFile = self.outputDir + "/sigmaPSF" + self.nr + ".txt"
        self.sim["PSF/AnalyticNonGaussian/Sigma/FromFile"] = psfInputFile
        myFile       = open(psfInputFile, "w")
        for sec in range(1000):
            myFile.write("{} {}\n".format(sec, 5 + sec / 25 ))
        myFile.close()







    def setEffectWithChargeDiffusion(self):
        # Sets the effects for the third run. Only one exposure is needed for this test.
        
        self.sim["ObservingParameters/NumExposures"]                = 1
        self.sim["PSF/AnalyticNonGaussian/ChargeDiffusionStrength"] = 10
        self.sim["PSF/AnalyticNonGaussian/IncludeChargeDiffusion"]  = "yes"
   






        
    def runSimulation(self):
        # Run the simulation when PSF sigma is constant.
        self.simFile1    = self.sim.run(removeOutputFile = True)

        # Run the simulation when PSF is given from file.
        self.setEffectForPSFFromFile()
        self.simFile2 = self.sim.run(removeOutputFile = True)

        # Run the simulation with charge diffusion.
        self.setEffectWithChargeDiffusion()
        self.simFile3 = self.sim.run(removeOutputFile = True)

        

        

        

    def compare(self):
        # Tests the validity of the simulations. The tests passes if all the sub-tests pass.
        
        simFile1    = self.simFile1
        simFile2    = self.simFile2
        simFile3    = self.simFile3
        
        condition1  = self.comparePositionCCD(simFile1)
        condition2  = self.comparePSFFromFile(simFile2)
        condition3  = self.compareDiffusion(simFile2, simFile3)

        return condition1 and condition2 and condition3
        







    def comparePositionCCD(self, simFile):
        # Tests that stars that are farther away from the optical axis are more spread out.

        xPos        = self.xPosition1
        yPos        = self.yPosition1
        image       = simFile.getImage(0)
        
        stars       = [np.array([[image[y+i-200][x+j-200] for i in range(400)] for j in range(400)]) for x, y in zip(xPos, yPos) ]

        xSpread     = [[star[i][i] for i in range(400)] for star in stars]
        ySpread     = [[star[i][399-i] for i in range(400)] for star in stars]

        sigma = self.getSigmas(stars)

        return self.doesIncreas(sigma)









    


    def comparePSFFromFile(self, simFile):
        # Test that the spread of the star is comparable with the one expected from the input file.
        
        numExpos = self.sim["ObservingParameters/NumExposures"]
        images   = [simFile.getImage(n) for n in range(numExpos)]
        sigma    = self.getSigmas(images)

        return self.doesIncreas(sigma)


    
    def compareDiffusion(self, simFile2, simFile3):
        # check that stars are more spread out when charge diffusion is included.
        
        image2 = simFile2.getImage(0)
        image3 = simFile3.getImage(0)
        images  = [image2, image3]

        sigma = self.getSigmas(images)

        return self.doesIncreas(sigma) 
 







    
    
    def getSigmas(self, spread, debug=False):
        # Estimates the sigma from a function.
 
        spread  = [ (star - np.min(star)) / (np.sum(star) - len(star) * np.min(star)) for star in spread]
        sigma = [fitGaussian2D(star, star.max(), len(spread)/2, len(spread)/2, 1, 1)[3]**2 + fitGaussian2D(star, star.max(), len(spread)/2, len(spread)/2, 1, 1)[4]**2 for star in spread]
        return sigma
    def doesIncreas(self, list):
        # Checks if a list in increasing.
        
        return all(np.array(list[:-1]) < np.array(list[1:]))


        











if __name__ == "__main__":

    t = AnalyticNonGaussianPSF()
    print(t.run())

    
