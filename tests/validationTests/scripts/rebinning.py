import numpy as np
import referenceFrames as rf
from test import Test
from math import degrees



"""
This test checks the debinnng done by PlatoSIM. It compares the result we obtained from PlatoSIM
from the one obtained by debinning using numpy functionality.
"""


class Rebinning(Test):

    def setNr(self):

        self.nr = "010"

    def setAllEffects(self):

        super().setAllEffects()

        self.sim["ObservingParameters/NumExposures"] = 1
        self.sim["SubField/NumRows"]    = 1500
        self.sim["SubField/NumColumns"] = 1500

        self.sim["PSF/Model"] = "MappedGaussian"
        self.sim["ControlHDF5Content/WriteSubPixelImages"] = "yes"


        # Create a SkyMap with a star that would fall onto that point in the sub field.
        ra, dec = rf.pixelToSkyCoordinates(self.sim, "2", 1000, 1000)

        ra  = degrees(ra)
        dec = degrees(dec)


        starCatalogFilename = self.outputDir + "/starCatalog"+ self.nr + ".txt"
        myFile = open(starCatalogFilename, "w")
        myFile.write("# RA DEC Vmag starID\n")
        myFile.write("{0}  {1}  {2}  {3}\n".format(ra, dec, 16.5, 1))
        myFile.close()

        self.sim["ObservingParameters/StarCatalogFile"] = starCatalogFilename



    def compare(self):

        image          = self.simFile.getImage(0)
        afterRebinning = self.rebinnedImage()

        return np.max(np.abs(afterRebinning - image)) * 100 < 0.001



    def rebinnedImage(self):

        imageBefore = self.simFile.getSubPixelImage(0)
        dimC        = self.sim["SubField/NumRows"]
        dimR        = self.sim["SubField/NumColumns"]
        numSubPix   = self.sim["SubField/SubPixels"]

        imageAfter  = imageBefore.reshape(dimC, numSubPix, dimR, numSubPix)
        return imageAfter.sum(axis=3).sum(axis=1)




if __name__ == "__main__":
    t = Rebinning()
    print(t.run())
    
