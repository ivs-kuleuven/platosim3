import random
import referenceFrames   as rf
import numpy             as np
import scipy.constants   as constants

from test import Test
from math import degrees


"""
This test is designed to check the field distortion. The test generates 20 start and runs the simulator with and without field distortion. The test checks weather the fields distortion is independent of the angle and heavily correlated with the radial distance from the optical axis.
"""








class FieldDistortion(Test):

    def setNr(self):

        self.nr = "004"



    def setAllEffects(self):

        super().setAllEffects()

        # Configure the input file and take only one exposure
        self.sim["ObservingParameters/NumExposures"] = 1
        self.sim["SubField/NumRows"] = 1000
        self.sim["SubField/NumColumns"] = 1000


        # Pick 100 random point in the sub field

        nRows    = self.sim["SubField/NumRows"]
        nColumns = self.sim["SubField/NumColumns"]

        self.pRows = [random.randint(1, nRows-1) for i in range(100)]
        self.pColumns = [random.randint(1, nColumns-1) for i in range(100)]

        # Create a SkyMap with a 20 stars in the sub field.

        starCatalogFilename = self.outputDir + "/starCatalog" + self.nr + ".txt"
        ra, dec = zip(*[rf.pixelToSkyCoordinates(self.sim, "2", row, col) for row, col in zip(self.pRows, self.pColumns)])

        ra  = [degrees(ra) for ra in ra]
        dec = [degrees(dec) for dec in dec]


        myFile = open(starCatalogFilename, "w")
        myFile.write("# RA DEC Vmag starID\n")
        i = 1
        for ra, dec in zip(ra, dec):
            magn = round(random.uniform(-5., 50.), 3)
            myFile.write("{0}  {1}  {2}  {3}\n".format(ra, dec, magn, i))
            i += 1
        myFile.close()

        self.sim["ObservingParameters/StarCatalogFile"] = starCatalogFilename






    def runSimulation(self):

        # Run the simulation without field distortion
        self.sim["Camera/IncludeFieldDistortion"] = "no"
        output = self.sim.run(removeOutputFile = True)

        self.resultUndistorted = output.getStarCoordinates(0)


        # Run the simulation with field distortion
        self.sim["Camera/IncludeFieldDistortion"] = "yes"
        output = self.sim.run(removeOutputFile = True)

        self.resultDistorted = output.getStarCoordinates(0)






    def compare(self):


        idUndistorted     = self.resultUndistorted[0]
        xUndistorted      = self.resultUndistorted[3]
        yUndistorted      = self.resultUndistorted[4]

        idDistorted       = self.resultDistorted[0]
        xDistorted        = self.resultDistorted[3]
        yDistorted        = self.resultDistorted[4]

        # Select those starts that are in both fields
        commonId          = np.intersect1d(idDistorted, idUndistorted)

        selectDistorted   = np.isin(idDistorted, commonId)
        selectUndistorted = np.isin(idUndistorted, commonId)

        xDistorted        = xDistorted[selectDistorted]
        yDistorted        = yDistorted[selectDistorted]


        xUndistorted      = xUndistorted[selectUndistorted]
        yUndistorted      = yUndistorted[selectUndistorted]

        # As measure for the field distortion we take the distance between
        # the undistorted and the distorted star position.
        xDelta = xDistorted - xUndistorted
        yDelta = yDistorted - yUndistorted

        dDist    = np.sqrt(xDelta**2 + yDelta**2)
        mDist    = np.mean(dDist)

        # The distortion is compared versus two independent parameters:
        # 1. The radial position of the undistorted stars.
        # 2. The angle of the undistorted stars.
        rUndist  = np.sqrt(xUndistorted**2 + yUndistorted**2)
        mUndist  = np.mean(rUndist)

        tUndist  = np.arctan(xUndistorted / yUndistorted)
        mtUndist = np.mean(tUndist)

        # The Pearson correlation coefficients between the distortion and radial distance,
        # and between the distortion and the angular position is calculated.
        corrRadN = np.dot(dDist - mDist, rUndist - mUndist)
        corrRadD = np.sqrt( (dDist - mDist).dot(dDist - mDist) * (rUndist - mUndist).dot(rUndist - mUndist) )
        corrR    = corrRadN / corrRadD

        corrAngN = np.dot(dDist - mDist, tUndist - mtUndist)
        corrAngD = np.sqrt((dDist - mDist).dot(dDist - mDist) * (tUndist - mtUndist).dot(tUndist - mtUndist))
        corrA    = corrAngN / corrAngD

        # The test passes if the correlation between the distortion and the radial distance is > 0.70 (strong
        # correlation) and the correlation between the distortion and the angular position < 0.3 (low correlation)
        #print(corrR, corrA)
        return (abs(corrR) > 0.7) and (corrA < 0.3)








if __name__ == "__main__":
    t = FieldDistortion()
    print(t.run())
