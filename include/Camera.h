#ifndef CAMERA_H
#define CAMERA_H

#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <array>

#include "Logger.h"
#include "Units.h"
#include "Constants.h"
#include "ConfigurationParameters.h"
#include "PointSpreadFunction.h"
#include "Polynomial1D.h"
#include "Heartbeat.h"
#include "HDF5File.h"
#include "HDF5Writer.h"
#include "Telescope.h"
#include "Detector.h"
#include "Sky.h"


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


    protected:

        Telescope &telescope;
        Sky &sky;

        double plateScale;                    // [arcsec/micron]
        double focalLength;                   // [mm]
        double focalPlaneOrientation;         // [rad]
        double throughputBandwidth;           // FWHM of the throughput passband [nm]
        double throughputLambdaC;             // Central wavelength of the throughput passband [nm]

        void selectPsf(double raStar, double decStar);
        
        pair<double, double> skyToAngularFocalPlaneCoordinates(double raStar, double decStar);
        pair<double, double> angularFocalPlaneToSkyCoordinates(double xFPprime, double yFPprime);

        pair<double, double> angularToPlanarFocalPlaneCoordinates(double xFPrad, double yFPrad);
        pair<double, double> planarToAngularFocalPlaneCoordinates(double xFPmm, double yFPmm);

        pair<double, double> planarToDistortedFocalPlaneCoordinates(double xFPmm, double yFPmm);

        double getGnomonicRadialDistanceFromOpticalAxis(double xFPprime, double yFPprime);

        void setDistortionPolynomial(Polynomial1D &polynomial);


    private:

        double internalTime;
        string polynomialType;
        double polynomialDegree;
        vector<double> polynomialCoefficients;

        PointSpreadFunction *psf;
        Polynomial1D polynomial;

        map<unsigned int, map<double, array<double, 6>>> detectedStarInfo;  // detectedStarInfo[starID][startTime] contains the values (xFPmean, yFPmean, sumFlux, Ndetections} values

};



#endif
