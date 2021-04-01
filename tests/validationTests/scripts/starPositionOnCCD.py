import random
import platosim.referenceFrames as rf
import numpy as np

from math import radians, degrees
from test import Test



"""
This test checks the position of a star on the CCD. It generates a random point on the CCD and based on that point
als a star map with a star that should fall onto this point. The test check weather the output file generates a
star on that position on the CCD +-1 pixel.
"""










class StarPositionOnCCD(Test):

    def setNr(self):

        self.nr = "001"


    def setAllEffects(self):

        # Configure the input file and take only one exposure
        super().setAllEffects()
        self.sim["ObservingParameters/NumExposures"] = 1

        # Pick a random point in the sub field
        nRows    = self.sim["SubField/NumRows"]
        nColumns = self.sim["SubField/NumColumns"]

        self.pRows    = random.randint(1, nRows-1)
        self.pColumns = random.randint(1, nColumns-1)


        # Create a SkyMap with a star that would fall onto that point in the sub field.
        ra, dec = rf.pixelToSkyCoordinates(self.sim, "2", self.pRows, self.pColumns)

        ra  = degrees(ra)
        dec = degrees(dec)


        starCatalogFilename = self.outputDir + "/starCatalog"+ self.nr + ".txt"
        myFile = open(starCatalogFilename, "w")
        myFile.write("# RA DEC Vmag starID\n")
        myFile.write("{0}  {1}  {2}  {3}\n".format(ra, dec, 16.5, 1))
        myFile.close()

        self.sim["ObservingParameters/StarCatalogFile"] = starCatalogFilename


    def compare(self):

        # Make sure the position on of the star in pixel coordinates is the same as the random
        # point we picked on the sub field. (with an error of +/- 1 pixel because of rounding
        # errors)
        image     = self.simFile.getImage(0)
        col , row = np.where(image == image.max())
        succes = (row[0] <= self.pRows + 1 <= row[0] + 2) and (col[0] <= self.pColumns + 1 <= col[0] + 2)

        return succes


if __name__ == "__main__":
    t = StarPositionOnCCD()
    print(t.run())
