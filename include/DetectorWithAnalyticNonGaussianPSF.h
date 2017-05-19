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




class IntegralOfAnalyticPSF 
{
    public:

        IntegralOfAnalyticPSF(size_t s) : size(s), n(0.) {}
        IntegralOfAnalyticPSF& addPart(double, double, double, double, double = 0., double = 0., double = 0.);
        double operator()(unsigned, unsigned, bool = true);

    private:

        size_t size;                              // number of (sub)pixels in one dimension
        double n;                                 // normalization factor
        vector<valarray<double>> erfxr;           // evaluated error functions for x
        vector<valarray<double>> erfyr;           // evaluated error functions for y
        vector<valarray<complex<double>>> erfxc;  // evaluated complex error functions for x
        vector<valarray<complex<double>>> erfyc;  // evaluated complex error functions for y
};










class DetectorWithAnalyticNonGaussianPSF: public Detector 
{
    public:

        DetectorWithAnalyticNonGaussianPSF(ConfigurationParameters &configParam, HDF5File &hdf5File, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator);
        virtual ~DetectorWithAnalyticNonGaussianPSF();

        virtual double takeExposure(int exposureNr, double startTime, double exposureTime) override;

        void configure(ConfigurationParameters &configParam);

        virtual tuple<bool, double, double> addFlux(double xFP, double yFP, double flux) override;
        virtual void addFlux(double flux) override;

        void integrateAnalyticPSF(IntegralOfAnalyticPSF&, double, double, double, double, double = 1.);

    protected:

        virtual void reset();
        virtual void integrateLight(int exposureNr, double startTime, double exposureTime) override;
        virtual void applyFlatfield() override;
        virtual void generateFlatfieldMap();
        virtual bool isInPixelMap(double row, double column);

        double sigma;                       // Width of the analytic PSF, equal to sigma for a Gaussian PSF
        vector<vector<double>> params;      // Table of analytic PSF parameters

        arma::Mat<float> flatfieldMap;      // Pixel flatfield map

        double flatfieldNoiseAmplitude;     // Peak-to-peak noise amplitude

        bool includeFlatfield;              // Whether or not to include flat fielding        
        long flatfieldSeed;

    private:

};


#endif
