#ifndef CAMERA_H
#define CAMERA_H

#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <array>

#include "armadillo"

#include "ArrayOperations.h"
#include "ConfigurationParameters.h"
#include "Constants.h"
#include "Detector.h"
#include "HDF5File.h"
#include "HDF5Writer.h"
#include "Heartbeat.h"
#include "Logger.h"
#include "Polynomial1D.h"
#include "Sky.h"
#include "StringUtilities.h"
#include "Platform.h"
#include "Telescope.h"
#include "Units.h"


using namespace std;


class Detector;  // forward declaration


class Camera : public HDF5Writer
{
    public:

        Camera(ConfigurationParameters &configParam, HDF5File &hdf5File, Platform &platform, Telescope &telescope, Sky &sky);
        virtual ~Camera();

        virtual void configure(ConfigurationParameters &configParam);
        virtual void exposeDetector(Detector &detector, double startTime, double exposureTime);

        virtual void initHDF5Groups() override;
        virtual void flushOutput() override;

        pair<double, double> skyToFocalPlaneCoordinates(double raStar, double decStar, bool useInitialOrientation=false);
        pair<double, double> focalPlaneToSkyCoordinates(double xFP, double yFP, bool useInitialOrientation=false);

        pair<double, double> undistortedToDistortedFocalPlaneCoordinates(double xFPmm, double yFPmm);
        pair<double, double> distortedToUndistortedFocalPlaneCoordinates(double xFPdist, double yFPdist);

        double getGnomonicRadialDistanceFromOpticalAxis(double xFP, double yFP);

        set<unsigned int> getAllStarIDs();

        double getTotalSkyBackground();


    protected:

        int beginExposureNr;                 // sequential number of first exposure. useful for slurm parallellisation
        int numExposures;                    // Number of exposures

        Platform &platform;
        Telescope &telescope;
        Sky &sky;

        double plateScale;                    // [arcsec/micron]
        double focalLength;                   // [mm]
        double throughputBandwidth;           // FWHM of the throughput passband [nm]
        double throughputLambdaC;             // Central wavelength of the throughput passband [nm]

        void setDistortionPolynomial(Polynomial1D &polynomial, Polynomial1D &inversePolynomial);

        double internalTime;
        string polynomialType;
        double polynomialDegree;
        vector<double> polynomialCoefficients;
        vector<double> inversePolynomialCoefficients;

        Polynomial1D polynomial;
        Polynomial1D inversePolynomial;

        bool includeFieldDistortion;      // Wheter or not field distortion should be included

        double userGivenSkyBackground;    // User-set zodiacal + stellar sky background.                        [phot/pix/s]
                                          // If negative, computed by the Sky class
        double fluxOfV0Star;              // Photon flux of a V=0 (G2V) star                                    [phot/s/m^2/nm]

        double raSun;                     // Right ascension of the direction of the sun shield during the run  [rad]
        double decSun;                    // Declination of the direction of the sun shield during the run      [rad]


        // detectedStarInfo[startTime][starID] contains the values (xFPmean, yFPmean, rowPixMean, colPixmean, sumFlux, Ndetections)

        map<double, map<unsigned int, array<double, 6>>> detectedStarInfo;
        vector<double> skyBackgroundValues;
        double totalSkyBackground;			// Total sky background [photons / pixel / exposure]

    private:

        
};



#endif
