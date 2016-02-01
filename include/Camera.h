#ifndef CAMERA_H
#define CAMERA_H

#include <string>
#include <cmath>

#include "Logger.h"
#include "TimeTicker.h"
#include "HDF5File.h"
#include "HDF5Writer.h"
#include "Detector.h"
#include "ConfigurationParameters.h"


using namespace std;



class Camera : public TimeTicker, HDF5Writer
{
    public:

        Camera(ConfigurationParameters &configParam, HDF5File &hdf5file);
        ~Camera();

        void exposeSubField(Detector &detector);

    protected:


    private:

        double plateScale;             // [arcsec/micron]
        double focalPlaneOrientation;  // [degrees]
        double internalTime;           // [s]

        void configure(ConfigurationParameters &configParam);
        void selectPsf(double raStar, double decStar);
        pair<double, double> getFocalPlaneCoordinates(double raStar, double decStar, 
            double raOpticalAxis, double decOpticalAxis, double plateScale);

};



#endif
