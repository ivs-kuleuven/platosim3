#ifndef DETECTORWITHANALYTICNONGAUSSIANPSF_H
#define DETECTORWITHANALYTICNONGAUSSIANPSF_H

#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <functional>
#include <complex>
#include <string>

#include "armadillo"

#include "Constants.h"
#include "Units.h"
#include "Detector.h"
#include "Camera.h"
#include "Parameter.h"
#include "Photometry.h"


using namespace std;



class DetectorWithAnalyticNonGaussianPSF: public Detector 
{
    public:

        DetectorWithAnalyticNonGaussianPSF(ConfigurationParameters &configParam,
					   HDF5File &hdf5File,
					   Camera &camera,
					   TemperatureGenerator &feeTemperatureGenerator,
					   TemperatureGenerator &detectorTemperatureGenerator,
					   Photometry &photometry,
					   double readoutTimeBeforeNextExposure,
					   double readoutTimeDuringNextExposure);
        ~DetectorWithAnalyticNonGaussianPSF();

        double takeExposure(int exposureNr, double startTime, double exposureTime) override;

        void configure(ConfigurationParameters &configParam);
        void updateParameters(double time) override;
        bool addFluxToMap(arma::Mat<float>& map, double row0, double col0, double r, double p, double flux);

        tuple<bool, double, double> addFlux(double xFP, double yFP, double flux) override;
        void addFlux(double flux) override;
        tuple<bool, double, double> addExtendedGhost(double xFP, double yFP, double radius, double flux) override;

        void integrateAnalyticPSF(IntegralOfAnalyticSignalResponse&, double, double, double, double, double, int = 1);
        void makeHighResolutionPSF(arma::Mat<float> &highResMap, bool includeDiffusion, int Npixels, int Nsubpixels);
        void generateDiffusionKernel(double kernelWidth);
        void applyDiffusionKernelOnPSF(double subpixRow, double subpixColumn, double flux, arma::fmat& psf, int numberOfPsfSubpixelsPerPixel);
        void applyPhotometry(const unsigned int exposureNr);
        void flushOutput() override;

    protected:

        void integrateLight(int exposureNr, double startTicme, double exposureTime) override;
        void applyFlatfield() override;
        void generateFlatfieldMap();

        Parameter<double> *sigma;                         // Width of the analytic PSF, equal to sigma for a Gaussian PSF
        vector<vector<double>> params;                    // Table of analytic PSF parameters
        arma::Mat<float> flatfieldMap;                    // Pixel flatfield map
        unsigned int numExposures;                        // Number of exposures
        unsigned int beginExposureNr;                     // Exposure nr of the first exposure in the time series
        double cycleTime;                                 // Image cycle time (exposure + readout before next exposure starts)  [s]
        double chargeDiffusionStrength;	                  // Strength of the charge diffusion (width of Gaussian diffusion kernel) [pixels]
        bool includeChargeDiffusion;	                  // Whether or not to include charge diffusion
        double flatfieldNoiseRMS;                         // Peak-to-peak noise amplitude
        bool includeFlatfield;                            // Whether or not to include flat fielding        
        long flatfieldSeed;                               // Seed dedicated to generate a random flatfield map
        bool writeFlatfieldMap;                           // Whether or not to write the flatfield map to the HDF5 file
        bool writeHighResolutionPSF;                      // Wheter or not to write the high resosultion PSF to the HDF5 file
        bool writeDiffusedPSF;                            // Wheter or not to write the high resosultion PSF to the HDF5 file 
        arma::Mat<float> diffusionKernel;                 // Diffusion kernel image
        IntegralOfAnalyticSignalResponse signalResponse;  // Signal response
        double diffusionKernelWidth;                      // Width (sigma) of the Gaussian diffusion kernel [sub-pixels]
        int diffusionKernelImageSize;                     // Size of the diffusion kernel image [sub-pixels]

    private:

};


#endif
