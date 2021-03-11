import platosim.referenceFrames as rf
import numpy                    as np

from math import degrees
from test import Test




"""This test is designed to check the sky background. The test runs the simulation when no stars are exposed to it. From the output file of the simulation, the theoretically predicted sky background is then subtracted. The resulting image should be filled with values around zero.   """










class SkyBackGround(Test):

    def setNr(self):
        self.nr = "007"


    def setAllEffects(self):

        # Configure the input file and take only one exposure
        super().setAllEffects()
        self.sim["ObservingParameters/NumExposures"] = 1

        # Create a sky map
        nRow     = self.sim["SubField/NumRows"]
        nCol     = self.sim["SubField/NumColumns"]
        row, col = nRow / 2, nCol / 2


        ra, dec = rf.pixelToSkyCoordinates(self.sim, "2", row, col)

        ra  = degrees(ra)
        dec = degrees(dec)

        starCatalogFilename = self.outputDir + "/starCatalog"+ self.nr + ".txt"
        myFile = open(starCatalogFilename, "w")
        myFile.write("# RA DEC Vmag starID\n")
        myFile.write("{0}  {1}  {2}  {3}\n".format(ra, dec, 16.5, 1))
        myFile.close()

        self.sim["ObservingParameters/StarCatalogFile"] = starCatalogFilename

        # Make sure no star does not fall in the sub field.
        self.sim["CCD/Position"] = "3"



    def compare(self):
        # This checks that the output image is approximately constant and around the expected sky background.

        image = self.simFile.getImage(0)


        exposureTime = self.sim["ObservingParameters/CycleTime"] - self.sim.getReadoutTime()[0]
        theoBackground = self.sim["Sky/SkyBackground"] * exposureTime * self.sim["Telescope/TransmissionEfficiency/BOL"]

        difference = image - int(theoBackground)
        RMS = np.sqrt(np.sum(difference*difference) / len(difference))
        return RMS < 0.01


if __name__ == "__main__":

    t = SkyBackGround()
    print(t.run())
