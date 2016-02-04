#ifndef CAMERA_H
#define CAMERA_H

#include <string>
#include <cmath>

#include "Logger.h"
#include "Units.h"
#include "Constants.h"
#include "ConfigurationParameters.h"
#include "HDF5File.h"
#include "Heartbeat.h"
#include "HDF5File.h"
#include "HDF5Writer.h"
#include "Telescope.h"
#include "Detector.h"
#include "Sky.h"
//#include "Telescope.h"


using namespace std;



class Camera : public HDF5Writer
{
    public:

        Camera(ConfigurationParameters &configParam, HDF5File &hdf5File, Telescope &telescope, Sky &sky);
        ~Camera();

        virtual void configure(ConfigurationParameters &configParam);
        virtual void exposeDetector(Detector &detector);

#ifdef UNIT_TESTS
        pair<double, double> test_skyToFocalPlaneCoordinates(double raStar, double decStar) {return skyToFocalPlaneCoordinates(raStar, decStar);};
        pair<double, double> test_focalPlaneToSkyCoordinates(double x, double y) {return focalPlaneToSkyCoordinates(x, y);};

        pair<double, double> test_convertOffsetMillimetersToDegrees(double xFP, double yFP) {return convertOffsetMillimetersToDegrees(xFP, yFP);};
        double test_getGnomonicRadialCoordinate(double xDeg, double yDeg) {return getGnomonicRadialCoordinate(xDeg, yDeg);};
#endif

    protected:

        Telescope &telescope;
        Sky &sky;

        double plateScale;                    // [arcsec/mm]
        double focalPlaneOrientation;         // [rad]
        double internalTime;                  // [s]

        void selectPsf(double raStar, double decStar);
        pair<double, double> skyToFocalPlaneCoordinates(double raStar, double decStar);
        pair<double, double> focalPlaneToSkyCoordinates(double x, double y);

        pair<double, double> convertOffsetMillimetersToDegrees(double xFP, double yFP);
        double getGnomonicRadialCoordinate(double xDeg, double yDeg);


    private:

};



#endif
