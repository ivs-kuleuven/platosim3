#ifndef DETECTORWITHSYMMETRICALMAPPEDPSF_H
#define DETECTORWITHSYMMETRICALMAPPEDPSF_H

#include <string>
#include <cmath>
#include <random>
#include <functional>

#include "armadillo"

#include "DetectorWithMappedPSF.h"
#include "SymmetricalPointSpreadFunction.h"

using namespace std;



class DetectorWithSymmetricalMappedPSF : public DetectorWithMappedPSF 
{
    public:

        DetectorWithSymmetricalMappedPSF(ConfigurationParameters &configParam, HDF5File &hdf5File, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure);
        virtual ~DetectorWithSymmetricalMappedPSF();

        void configure(ConfigurationParameters &configParam);

    protected:

        SymmetricalPointSpreadFunction *psf;

        void setPsfForSubfield();

};

#endif
