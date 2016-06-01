#ifndef CAMERA_H
#define CAMERA_H

#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <array>

#include "ArrayOperations.h"
#include "ConfigurationParameters.h"
#include "Constants.h"
#include "Detector.h"
#include "HDF5File.h"
#include "HDF5Writer.h"
#include "Heartbeat.h"
#include "Logger.h"
#include "PointSpreadFunction.h"
#include "Polynomial1D.h"
#include "Sky.h"
#include "StringUtilities.h"
#include "Telescope.h"
#include "Units.h"


using namespace std;


class Detector;  // forward declaration


class Camera : public HDF5Writer
{
    public:

        Camera(ConfigurationParameters &configParam, HDF5File &hdf5File, Telescope &telescope, Sky &sky);
        virtual ~Camera();

        virtual void configure(ConfigurationParameters &configParam);
        virtual void exposeDetector(Detector &detector, double startTime, double exposureTime);

        virtual void initHDF5Groups() override;
        virtual void flushOutput() override;

        virtual arma::fmat getRebinnedPsfForPlanarFocalPlaneCoordinates(double xFPmm, double yFPmm, unsigned int targetSubPixels, double orientationAngle);

        pair<double, double> skyToAngularFocalPlaneCoordinates(double raStar, double decStar);
        pair<double, double> angularFocalPlaneToSkyCoordinates(double xFPprime, double yFPprime);

        pair<double, double> angularToPlanarFocalPlaneCoordinates(double xFPrad, double yFPrad);
        pair<double, double> planarToAngularFocalPlaneCoordinates(double xFPmm, double yFPmm);

        pair<double, double> planarToDistortedFocalPlaneCoordinates(double xFPmm, double yFPmm);
        pair<double, double> distortedToPlanarFocalPlaneCoordinates(double xFPdist, double yFPdist);

        double getGnomonicRadialDistanceFromOpticalAxis(double xFPprime, double yFPprime);

        set<unsigned int> getAllStarIDs();


    protected:

        Telescope &telescope;
        Sky &sky;

        double plateScale;                    // [arcsec/micron]
        double focalLength;                   // [mm]
        double focalPlaneOrientation;         // [rad]
        double throughputBandwidth;           // FWHM of the throughput passband [nm]
        double throughputLambdaC;             // Central wavelength of the throughput passband [nm]

        void setDistortionPolynomial(Polynomial1D &polynomial, Polynomial1D &inversePolynomial);


    private:

        double internalTime;
        string polynomialType;
        double polynomialDegree;
        vector<double> polynomialCoefficients;
        vector<double> inversePolynomialCoefficients;

        PointSpreadFunction *psf;
        Polynomial1D polynomial;
        Polynomial1D inversePolynomial;

        bool includeFieldDistortion;          // Wheter or not field distortion should be included

        double userGivenSkyBackground;        // User-set zodiacal + stellar sky background. [phot/pix/s]
                                              // If negative, computed by the Sky class

        double fluxOfV0Star;                  // Photon flux of a V=0 (G2V) star [phot/s/m^2/nm]

        // detectedStarInfo[startTime][starID] contains the values 
        //    (xFPmean, yFPmean, rowPixMean, colPixmean, sumFlux, Ndetections)

        map<double, map<unsigned int, array<double, 6>>> detectedStarInfo;
        vector<double> skyBackgroundValues;

};



#endif
