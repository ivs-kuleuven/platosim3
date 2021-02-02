import numpy           as np
import referenceFrames as rf
from test       import Test
from math       import radians
from validation import equatorial2galactic, galactic2equatorial, aberration
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


        dim = 6
        self.sim["SubField/NumRows"]         = dim
        self.sim["SubField/NumColumns"]      = dim
        self.sim["SubField/ZeroPointRow"]    = 100 - dim // 2
        self.sim["SubField/ZeroPointColumn"] = 100 - dim // 2




        self.sim["ObservingParameters/NumExposures"]   = 1
        self.sim["Camera/IncludeAberrationCorrection"] = "yes"
        self.sim["Camera/AberrationCorrection/Type"]   = "absolute"


    def runSimulation(self):


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

        for exp in range(0, numExposures, deltaNumExposures):
            self.sim["ObservingParameters/BeginExposureNr"] = exp
            self.simFile = self.sim.run(removeOutputFile=True)
            ([x], [y])   = self.simFile.getStarCoordinates(exp)[3:5]
            ra, dec      = rf.focalPlaneToSkyCoordinates(x, y, raPlatform,
                decPlatform, solarPanelOrientation, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength)
            outputRa.append(ra)
            outputDec.append(dec)

        self.outputRa  = np.rad2deg(np.array(outputRa))
        self.outputDec = np.rad2deg(np.array(outputDec))


    def compare(self):


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
        radialTheo = [x**2 + y**2 for x, y in zip(dxTheo, dyTheo)]

        # Check that the observed displacements fall within the theoretical
        #predicted range.
        rTheoMin   = math.sqrt(min(radialTheo))
        rObseMin   = math.sqrt(min(radialDistance))

        rTheoMax   = math.sqrt(max(radialTheo))
        rObseMax   = math.sqrt(max(radialDistance))

        condition1 = (rTheoMin - rObseMin) / rObseMin < 0.01
        condition2 = (rObseMax - rTheoMax) / rObseMax < 0.01

        return ( condition1 and condition2 )




#t3 = AbsoluteAberration()
#t3.run()
