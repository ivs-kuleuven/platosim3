#ifndef CAMERA_H
#define CAMERA_H

#include <string>
#include <cmath>
#include <algorithm>

#include "Logger.h"
#include "Units.h"
#include "Constants.h"
#include "ConfigurationParameters.h"
#include "PointSpreadFunction.h"
#include "HDF5File.h"
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
        ~Camera();

        virtual void configure(ConfigurationParameters &configParam);
        virtual void exposeDetector(Detector &detector, double startTime, double exposureTime);

    protected:

        Telescope &telescope;
        Sky &sky;

        double plateScale;                    // [arcsec/mm]
        double focalPlaneOrientation;         // [rad]
        double throughputBandwidth;           // FWHM of the throughput passband [nm]
        double throughputLambdaC;             // Central wavelength of the throughput passband [nm]

        void selectPsf(double raStar, double decStar);
        pair<double, double> skyToFocalPlaneCoordinates(double raStar, double decStar);
        pair<double, double> skyToNormalizedFocalPlaneCoordinates(double raStar, double decStar);
        pair<double, double> focalPlaneToSkyCoordinates(double x, double y);

        double getGnomonicRadialDistanceFromOpticalAxis(double xFPprime, double yFPprime);
        double getGnomonicRadialDistanceFromOpticalAxisNormalized(double xFPprime, double yFPprime);

    private:

        double internalTime;

        PointSpreadFunction *psf;

};



#endif
