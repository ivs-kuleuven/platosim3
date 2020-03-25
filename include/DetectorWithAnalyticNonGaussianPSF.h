#ifndef DETECTORWITHANALYTICNONGAUSSIANPSF_H
#define DETECTORWITHANALYTICNONGAUSSIANPSF_H

#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <functional>
#include <complex>

#include "armadillo"

#include "Constants.h"
#include "Units.h"
#include "Detector.h"
#include "Parameter.h"


using namespace std;



class DetectorWithAnalyticNonGaussianPSF: public Detector 
{
    public:

        DetectorWithAnalyticNonGaussianPSF(ConfigurationParameters &configParam, HDF5File &hdf5File, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure);
        virtual ~DetectorWithAnalyticNonGaussianPSF();

        virtual double takeExposure(int exposureNr, double startTime, double exposureTime) override;

        void configure(ConfigurationParameters &configParam);
        virtual void updateParameters(double time) override;

        virtual tuple<bool, double, double> addFlux(double xFP, double yFP, double flux, int subsubfieldx, int subsubfieldy) override;  //%% Added subsubfields for spectral dependency
        virtual void addFlux(double flux, int subsubfieldx, int subsubfieldy) override;  //%% Added subsubfields for spectral dependency

        void integrateAnalyticPSF(IntegralOfAnalyticSignalResponse&, double, double, double, double, int = 1);

        void makeHighResolutionPSF(arma::Mat<float> &highResMap, int Npixels, int Nsubpixels);
        virtual void flushOutput() override;

    protected:

        virtual void reset();
        virtual void integrateLight(int exposureNr, double startTime, double exposureTime, int subsubfieldx, int subsubfieldy) override;  //%% Added subsubfields for spectral dependency
        virtual void applyFlatfield(int subsubfieldx, int subsubfieldy) override;  //%% Added subsubfields for spectral dependency
        virtual void generateFlatfieldMap();

        Parameter<double> *sigma;           // Width of the analytic PSF, equal to sigma for a Gaussian PSF
        vector<vector<double>> params;      // Table of analytic PSF parameters

        arma::Mat<float> flatfieldMap;      // Pixel flatfield map

        double chargeDiffusionStrength;			// Strength of the charge diffusion (width of the Gaussian diffusion kernel) [pixels]
        bool includeChargeDiffusion;				// Whether or not to include charge diffusion

        double flatfieldNoiseRMS;     // Peak-to-peak noise amplitude

        int wave_bins;  //%%  Number of wavelength bins to be processed, added for spectral dependency

        bool includeFlatfield;              // Whether or not to include flat fielding        
        long flatfieldSeed;

    private:

};


#endif
