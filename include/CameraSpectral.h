#ifndef CAMERASPECTRAL_H
#define CAMERASPECTRAL_H

#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <array>
#include <tuple>

#include "armadillo"

#include "Camera.h"
#include "SpectralDependenceUtility.h"


using namespace std;


class Detector;  // forward declaration

typedef map<unsigned int, array<double, 6>>::iterator starInfoIterator;


class CameraSpectral : public Camera
{
    public:

        CameraSpectral(ConfigurationParameters &configParam, HDF5File &hdf5File, Platform &platform, Telescope &telescope, Sky &sky);
        virtual ~CameraSpectral(){};

        virtual void exposeDetector(Detector &detector, double startTime, double exposureTime, double readoutTimeBeforeNextExposure);

    protected:

    private:

        double referenceWavelength;
        int binnumber;
        double lowerWavelength;
        double binwidth;
        vector<double> transmissionEfficiencySpectral;

        
};



#endif
