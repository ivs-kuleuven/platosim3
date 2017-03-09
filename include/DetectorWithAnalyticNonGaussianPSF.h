#ifndef DETECTORWITHANALYTICNONGAUSSIANPSF_H
#define DETECTORWITHANALYTICNONGAUSSIANPSF_H

#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <functional>
#include <complex>
#include <valarray>

#include "armadillo"
#include "Faddeeva.hh"

#include "Constants.h"
#include "Units.h"
#include "Detector.h"

using namespace std;



class DetectorWithAnalyticNonGaussianPSF: public Detector 
{
    public:

        DetectorWithAnalyticNonGaussianPSF(ConfigurationParameters &configParam, HDF5File &hdf5File, Camera &camera);
        virtual ~DetectorWithAnalyticNonGaussianPSF();

        virtual double takeExposure(int exposureNr, double startTime, double exposureTime) override;

        void configure(ConfigurationParameters &configParam);

        virtual tuple<bool, double, double> addFlux(double xFP, double yFP, double flux) override;
        virtual void addFlux(double flux) override;

    protected:

        virtual void reset();
        virtual void integrateLight(int exposureNr, double startTime, double exposureTime) override;
        virtual void applyFlatfield() override;
        virtual void generateFlatfieldMap();
        virtual bool isInPixelMap(double row, double column);

        arma::Mat<float> flatfieldMap;      // Pixel flatfield map

        double flatfieldNoiseAmplitude;     // Peak-to-peak noise amplitude

        bool includeFlatfield;              // Whether or not to include flat fielding        
        long flatfieldSeed;

    private:

};


#endif
