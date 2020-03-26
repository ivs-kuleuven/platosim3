#ifndef DETECTORWITHMAPPEDPSF_H
#define DETECTORWITHMAPPEDPSF_H

#include <string>
#include <cmath>
#include <random>
#include <functional>

#include "armadillo"

#include "Detector.h"
#include "PointSpreadFunction.h"

using namespace std;



class DetectorWithMappedPSF : public Detector 
{
    public:

        DetectorWithMappedPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure);
        virtual ~DetectorWithMappedPSF(){};

        virtual double takeExposure(int exposureNr, double startTime, double exposureTime) override;


        void configure(ConfigurationParameters &configParam);

        virtual tuple<bool, double, double> addFlux(double xFP, double yFP, double flux, int subsubfieldx, int subsubfieldy) override;  //%% Added subsubfields for spectral dependency
        virtual void addFlux(double flux, int subsubfieldx, int subsubfieldy) override;  //%% Added subsubfields for spectral dependency

//@        virtual void configure(ConfigurationParameters &configParam){};

    protected:

        virtual void initHDF5Groups() override;

        virtual void reset();
        
        virtual void integrateLight(int exposureNr, double startTime, double exposureTime, int susubfieldx, int subsubfieldy) override;  //%% Added subsubfields for spectral dependency
        virtual bool isInSubPixelMap(double row, double column);
        virtual void applyDiffusionKernel(double row, double column, double flux);
        virtual void applyFlatfield(int subsubfieldx, int subsubfieldy) override;  //%% Added subsubfields for spectral dependency

        virtual void generateFlatfieldMap();
        virtual void generateDiffusionKernel(double kernelWidth);
        virtual void rebin();
        void writeSubPixelMapToHDF5(int exposureNr);

        void setPsfForSubfield(int subsubfieldx, int subsubfieldy);  //%% Added subsubfields for spectral dependency
        virtual void convolveWithPsf(int binnumber);  //%% Added wavelength bin for spectral dependency

        arma::Mat<float> subPixelMap;           // Sub-pixel map, incl. edge pixels
        arma::Mat<float> psfMap;                // The PSF map that will be used for convolving

        vector<arma::Mat<float>> psfVector;  //%% Added for spectral dependency, vector of all psfMaps to be used

        arma::Mat<float> flatfieldMap;          // Intra-pixel flatfield map

        double chargeDiffusionStrength;			// Strength of the charge diffusion (width of the Gaussian diffusion kernel) [pixels]
        bool includeChargeDiffusion;			// Whether or not to include charge diffusion
        bool includeJitterSmoothing;            // Whether or not to include jitter smoothing

        double flatfieldNoiseRMS;               // Peak-to-peak noise amplitude

	int wave_bins;  //%%  Number of wavelength bins to be processed, added for spectral dependency

        bool includeFlatfield;                  // Whether or not to include flat fielding
        bool writeFlatfieldMap;                 // Whether or not to write the flatfield map to the HDF5 file
        bool writeSubPixelImagesToHDF5;         // Write subpixel maps to HDF5 as well
        bool includeConvolution;                // Whether or not to convolve the subPixelMap with the PSF

        arma::Mat<float> diffusionKernel;                    // Diffusion kernel image
        IntegralOfAnalyticSignalResponse signalResponse;	 // Signal response
        double diffusionKernelWidth;				         // Width (sigma) of the Gaussian diffusion kernel [sub-pixels]
        int diffusionKernelImageSize;             		     // Size of the diffusion kernel image [sub-pixels]

        unsigned int numRowsSubPixelMap;        // Nr of subpixel rows in the subfield incl. edge pixels (= size in the y-direction) [subpixels]
        unsigned int numColumnsSubPixelMap;     // Nr of subpixel columns in the subfield incl. edge pixels (= size in the x-direction = readout direction) [subpixels]
        unsigned int numSubPixelsPerPixel;      // Nr of sub-pixels per pixel
        
        long flatfieldSeed;

        Convolver convolver;
};

#endif
