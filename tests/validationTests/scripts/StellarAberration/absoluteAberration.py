import numpy           as np
import platosim.referenceFrames as rf
from test       import Test
from math       import radians
from platosim.validation import equatorial2galactic, galactic2equatorial, aberration
import math
import matplotlib.pyplot as plt



class AbsoluteAberration(Test):

    def setNr(self):
        self.nr = "003.1"

    def setAllEffects(self):
        super().setAllEffects()

        # Define star catalog file
        starFileName = self.outputDir + "/starCatalog" + self.nr + ".txt"
        pos  = np.array([100])
        mag  = np.array([12.5])
        self.sim.createStarCatalogFileFromPixelCoordinates(pos, pos, mag, np.array([1]), starFileName)
        self.sim["ObservingParameters/StarCatalogFile"] = starFileName


        dim = 4510
        self.sim["SubField/NumRows"]         = dim
        self.sim["SubField/NumColumns"]      = dim
        self.sim["SubField/ZeroPointRow"]    = 0
        self.sim["SubField/ZeroPointColumn"] = 0


        self.sim["ObservingParameters/NumExposures"]   = 1
        self.sim["Camera/IncludeAberrationCorrection"] = "yes"
        self.sim["Camera/AberrationCorrection/Type"]   = "absolute"

        self.startTime = 9442527.53198146
        self.velocity  = [-0.4283823880272332, 0.9035748636371567, 0.006402767462494304]
        self.sim["Camera/AberrationCorrection/StartTime"] = self.startTime



    def runForYear(self):
        # We run the simulation multiple times over one year with 5 days
        #between two consecutive exposures.
        raPlatform            = radians(self.sim["ObservingParameters/RApointing"])
        decPlatform           = radians(self.sim["ObservingParameters/DecPointing"])
        solarPanelOrientation = radians(self.sim["Platform/SolarPanelOrientation"])
        tiltAngle             = radians(self.sim["Telescope/TiltAngle"])
        azimuthAngle          = radians(self.sim["Telescope/AzimuthAngle"])
        focalPlaneAngle       = radians(self.sim["Camera/FocalPlaneOrientation/ConstantValue"])
        focalLength           = self.sim["Camera/FocalLength/ConstantValue"] * 1000
        numExposures          = 365 * 24 * 60 * 60 // 25
        deltaNumExposures     = 24 * 60 * 60 // 25 * 5

        outputRa  = []
        outputDec = []
        startTime0 = self.startTime
        
        for exp in range(0, numExposures + 1, deltaNumExposures):
            self.sim["ObservingParameters/BeginExposureNr"] = exp
            self.sim["Camera/AberrationCorrection/StartTime"] = startTime0 + exp*25
            self.simFile = self.sim.run(removeOutputFile=True)
            ([x], [y])   = self.simFile.getStarCoordinates(exp)[3:5]
            ra, dec      = rf.focalPlaneToSkyCoordinates(x, y, raPlatform,
                decPlatform, solarPanelOrientation, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength)
            outputRa.append(ra)
            outputDec.append(dec)

        self.outputRa  = np.rad2deg(np.array(outputRa))
        self.outputDec = np.rad2deg(np.array(outputDec))



    def runForSecondTest(self):
        lamb = np.arctan2(self.velocity[1], self.velocity[0])
        beta = np.arcsin(self.velocity[2])
        ra, dec = np.rad2deg(rf.ecliptic2equatorial(lamb, beta))
        print("ra: ", ra, "dec: ",  dec)
        self.sim["ObservingParameters/RApointing"] = ra
        self.sim["ObservingParameters/DecPointing"] = dec
        self.sim["CCD/Position"] = "1"

        # Set the star catalog
        starFileName = self.outputDir + "/starCatalog" + self.nr + ".txt"
        myFile = open(starFileName, "w")
        myFile.write("# RA DEC Vmag starID\n")
        myFile.write("{0}  {1}  {2}  {3}\n".format(ra, dec, 16.5, 1))
        myFile.close()

        #print(rf.skyToFocalPlaneCoordinates(ra, dec, 0, 0, 0, 0, 0, 0, 247.52))
        


        
    def runSimulation(self):

        # Visually check that the abberation is in line with a theoretical estimate (with constant velocity 30km/s)
        self.runForYear()
        self.compareToTheoreticalEstimate()

        #self.runForSecondTest()

        


    def compareToTheoreticalEstimate(self):

        # The theoretical predicted aberration is calculated over the span of one year. We check whether the observed displacements over
        # one year fall withing the predicted range.
            
        ([inputRa], [inputDec]) = (self.simFile.getStarCatalog()[1:3])
        lonStar, latStar = equatorial2galactic(inputRa, inputDec)
            
        # Displacement that is measured over the span of a year.
        dxMeasurement  = (self.outputRa - inputRa) * 3600 * np.cos(np.deg2rad(self.outputDec))
        dyMeasurement  = (self.outputDec - inputDec) * 3600
        radialDistance = [x0**2 + y0**2 for x0, y0 in zip(dxMeasurement, dyMeasurement)]
            

        # Theoretically predicted aberration over one year
        theoAber = [aberration(lonStar, latStar, lSun) for lSun in np.linspace(0, 2*np.pi, 360)]
        tAberGal = [galactic2equatorial(lamStar, betaStar) for lamStar, betaStar in theoAber]
        
        raTheo, decTheo = zip(*tAberGal)
        
        dxTheo     = (raTheo - inputRa) * 3600 * np.cos(np.deg2rad(decTheo))
        dyTheo     = (decTheo - inputDec ) * 3600

        
        plt.plot(dxMeasurement, dyMeasurement, 'r')
        plt.plot(dxTheo, dyTheo)
        plt.savefig(self.ioPath + "/test" + self.nr +  "/aberrationOverYear.png")
        


    def compare(self):

        #result1 = self.NoAberration()
        return True



if __name__ == "__main__":
    t = AbsoluteAberration()
    print(t.run())
