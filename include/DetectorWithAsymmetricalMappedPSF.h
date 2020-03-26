#ifndef DETECTORWITHASYMMETRICALMAPPEDPSF_H
#define DETECTORWITHASYMMETRICALMAPPEDPSF_H

#include <string>
#include <cmath>
#include <random>
#include <functional>

#include "armadillo"

#include "DetectorWithMappedPSF.h"
#include "AsymmetricalPointSpreadFunction.h"

using namespace std;



class DetectorWithAsymmetricalMappedPSF : public DetectorWithMappedPSF 
{
    public:

        DetectorWithAsymmetricalMappedPSF(ConfigurationParameters &configParam, HDF5File &hdf5File, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure);
        virtual ~DetectorWithAsymmetricalMappedPSF();

        void configure(ConfigurationParameters &configParam);

    protected:

        AsymmetricalPointSpreadFunction *psf;

        void setPsfForSubfield(int subsubfieldx, int subsubfieldy);

};

#endif
