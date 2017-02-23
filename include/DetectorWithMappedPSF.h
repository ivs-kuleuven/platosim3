#ifndef DETECTORWITHMAPPEDPSF_H
#define DETECTORWITHMAPPEDPSF_H

#include <string>
#include <cmath>
#include <random>
#include <functional>

#include "armadillo"

#include "Detector.h"

using namespace std;



class DetectorWithMappedPSF: public Detector 
{
    public:

        DetectorWithMappedPSF(ConfigurationParameters &configParam, HDF5File &hdf5File, Camera &camera);
        virtual ~DetectorWithMappedPSF();

        virtual double takeExposure(double startTime, double exposureTime) override;

        void configure(ConfigurationParameters &configParam);

        virtual tuple<bool, double, double> addFlux(double xFP, double yFP, double flux) override;
        virtual void addFlux(double flux) override;

    protected:

        virtual void reset();
        virtual void initHDF5Groups() override;
        virtual void integrateLight(double startTime, double exposureTime) override;
        virtual bool isInSubPixelMap(double row, double column);
        virtual void applyFlatfield() override;
        virtual void generateFlatfieldMap();
        virtual void rebin();
        void writeSubPixelMapToHDF5();

        void setPsfForSubfield();
        virtual void convolveWithPsf();

        PointSpreadFunction *psf;

        arma::Mat<float> subPixelMap;           // Sub-pixel map, incl. edge pixels
        arma::Mat<float> psfMap;                // The PSF map that will be used for convolving
        arma::Mat<float> flatfieldMap;          // Intra-pixel flatfield map

        double flatfieldNoiseAmplitude;         // Peak-to-peak noise amplitude

        bool includeFlatfield;                  // Whether or not to include flat fielding
        bool writeSubPixelImagesToHDF5;         // Write subpixel maps to HDF5 as well
        bool includeConvolution;                // Whether or not to convolve the subPixelMap with the PSF


        unsigned int numRowsSubPixelMap;        // Nr of subpixel rows in the subfield incl. edge pixels (= size in the y-direction) [subpixels]
        unsigned int numColumnsSubPixelMap;     // Nr of subpixel columns in the subfield incl. edge pixels (= size in the x-direction = readout direction) [subpixels]
        unsigned int numSubPixelsPerPixel;      // Nr of sub-pixels per pixel
        
        long flatfieldSeed;



    private:

        Convolver convolver;
        int subPixelImageNr;

};


#endif
