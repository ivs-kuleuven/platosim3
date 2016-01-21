
#ifndef DETECTOR_H
#define DETECTOR_H

#include <string>
#include <cmath>
#include <random>
#include "logger.h"
#include "hdf5file.h"
#include "hdf5writer.h"


using namespace std;



class Detector : public HDF5Writer
{
    public:

        Detector(HDF5File &hdf5File);
        virtual ~Detector();

        virtual void takeExposure(double startTime, double exposureTime);

    protected:

        virtual void reset();
        virtual void generateFlatfieldMap();

        virtual void integrateLight(double startTime, double exposureTime);
        virtual bool isInSubPixelMap(double row, double column);
        virtual void addFlux(double xCoords, double yCoords, double flux);
        virtual void addFlux(double flux);
        virtual void applyFlatfield();
        virtual void rebin();
        
        virtual void readOut(double exposureTime);
        virtual void applyQuantumEfficiency();
    	virtual void addPhotonNoise();
    	virtual void applyFullWellSaturation();
    	virtual void applyCTE();
    	virtual void applyOpenShutterSmearing();
    	virtual void addReadoutNoise();
    	virtual void applyGain();
    	virtual void addElectronicOffset();	     
    	virtual void applyDigitalSaturation();
    
        virtual void initHDF5Groups() override;
        void writePixelMapToHDF5();


        arma::Mat<float> pixelMap;               // Pixel map, excl. edge pixels
        arma::Mat<float> subPixelMap;            // Sub-pixel map, incl. edge pixels
        arma::Mat<float> smearingMap;            // Smearing map (i.e. over-scan strip)
        arma::Mat<float> biasMap;                // Bias map (i.e. pre-scan strip)
        arma::Mat<float> cteMap;                 // CTE map
        arma::Mat<float> flatfieldMap;           // Flatfield map

        unsigned int numRows;                    // Nr of rows of the detector (= size in y-direction) [pixels]
    	unsigned int numColumns;                 // Nr of columns of the detector (= size in x-direction = readout direction) [pixels]
        unsigned int numRowsSubField;            // Nr of rows in the subfield excl. edge pixels (= size the y-direction) [pixels]
        unsigned int numColumnsSubField;         // Nr of columns in the subfield excl. edge pixels (= size in the x-direction = readout direction) [pixels]
        unsigned int numRowsSubPixelMap;         // Nr of subpixel rows in the subfield incl. edge pixels (= size in the y-direction) [subpixels]
        unsigned int numColumnsSubPixelMap;      // Nr of subpixel columns in the subfield incl. edge pixels (= size in the x-direction = readout direction) [subpixels]
        unsigned int numRowsSmearingMap;         // Nr of rows in the smearing overscan strip [pixels]
        unsigned int numRowsBiasMap;             // Nr of rows in the bias prescan strip [pixels]

    	double originOffsetRow;                  // Offset of the detector origin from the centre of the optical plane in the row direction [mm]
    	double originOffsetColumn;               // Offset of the detector origin from the centre of the optical plane in the column direction [mm]
        unsigned int subFieldZeroPointRow;       // Position of the subfield zeropoint w.r.t. the complete detector in the row direction [pixels]
        unsigned int subFieldZeroPointColumn;    // Position of the subfield zeropoint w.r.t. the complete detector in the column direction [pixels]
    	double orientationAngle;                 // Orientation angle of the detector w.r.t. the orientation of the focal plane, measured counterclockwise [degrees]
 
    	double pixelSize;	                     // Pixel size [microns]
        unsigned int numSubPixelsPerPixel;	     // Nr of sub-pixels per pixel
    	unsigned int numEdgePixels;              // Nr of pixels to extend the subfield on each side, to account for the edge effect


    	double flatfieldNoiseAmplitude;          // Peak-to-peak noise amplitude

    	double quantumEfficiency;	             // Quantum efficiency (in [0,1])
    	double readoutTime;                      // Readout time [s]
    	double chargeTransferTime;	             // Charge transfer time [s]
        double meanCte;                          // Mean charge-transfer efficiency
        double readoutNoise;                     // Mean readout noise [electrons]
        double gain;                             // Detector gain [electrons / ADU]
        unsigned long fullWellSaturationLimit;   // Full-well saturation limit [electrons/pixel]
        unsigned int electronicOffset;           // Bias or electronic offset [ADU]
        unsigned long digitalSaturationLimit;    // Digital saturation limit [ADU / pixel]

    	bool includePhotonNoise;                 // Whether or not to include photon noise

        double internalTime;

    	double flatfieldSeed;
    	double readoutNoiseSeed;
    	double photonNoiseSeed;
    	double cteMapSeedRow;
    	double cteMapSeedColumn;

    	mt19937 photonNoiseGenerator;
    	mt19937 readoutNoiseGenerator;

    	poisson_distribution<int> photonNoiseDistribution;
    	normal_distribution<double> readoutNoiseDistribution;



    private:

        int imageNr;

};



#endif
