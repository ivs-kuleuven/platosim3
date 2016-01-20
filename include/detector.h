
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

        virtual void initFlatfieldMap();

        virtual void reset();
    
        virtual void integrateLight(double startTime, double exposureTime);
        virtual void addFlux(double xCoords, double yCoords, double flux);
        virtual bool isInSubPixelMap(double row, double column);
        virtual void addFlux(double flux);
        virtual void applyFlatfield();
        virtual void rebin();
        
        virtual void readOut(double exposureTime);
        virtual void applyQuantumEfficiency();
    	virtual void addPhotonNoise();
    	virtual void applyFullWellSaturation();
    	virtual void applyCte();
    	virtual void applyOpenShutterSmearing();
    	virtual void addReadoutNoise();
    	virtual void applyGain();
    	virtual void addElectronicOffset();	     // Bias
    	virtual void applyDigitalSaturation();
    
        virtual void initHDF5Groups() override;

        void writePixelMapToHDF5();


        arma::Mat<float> pixelMap;               // Pixel map, excl. edge pixels

        unsigned int numRows;                    // Number of rows of the detector (i.e. dimension in the y-direction) [pixels]
    	unsigned int numColumns;                 // Number of columns of the detector (i.e. dimension in the x-direction = readout direction) [pixels]
    	double originOffsetRow;                  // Offset of the detector origin from the centre of the optical plane in the row direction [mm]
    	double originOffsetColumn;               // Offset of the detector origin from the centre of the optical plane in the column direction [mm]
    	double orientationAngle;                 // Orientation angle of the detector w.r.t. the orientation of the focal plane, measured counterclockwise [degrees]
 
    	unsigned int numRowsSmearingMap;         // Number of rows in the smearing over-scan strip [pixels]
    	unsigned int numRowsBiasMap;	         // Number of rows in the bias pre-scan strip [pixels]

    	double pixelSize;	                     // Pixel size [microns]
    

    	unsigned int subFieldZeroPointRow;	     // Position of the sub-field zeropoint w.r.t. the complete detector in the row direction [pixels]
    	unsigned int subFieldZeroPointColumn;    // Position of the sub-field zeropoint w.r.t. the complete detector in the column direction [pixels]
    	unsigned int numRowsSubField;	         // Number of rows in the sub-field at pixel level and excl. edge pixels (i.e. dimension in the y-direction) [pixels]
    	unsigned int numColumnsSubField;	     // Number of columns in the sub-field at pixel level and excl. edge pixels  (i.e. dimension in the x-direction = readout direction) [pixels]
    	unsigned int numSubPixelsPerPixel;	     // Number of sub-pixels per pixel
    	unsigned int numEdgePixels;              // Number of pixels to extend the sub-field on each side, to account for the edge effect

    	arma::Mat<float> subPixelMap;	         // Sub-pixel map, incl. edge pixels

    	unsigned int numRowsSubPixelMap;	     // Number of rows in the sub-field at sub-pixel level and incl. edge pixels (i.e. dimensions in the y-direction) [sub-pixels]
    	unsigned int numColumnsSubPixelMap;	     // Number of columns in the sub-field at sub-pixel level and incl. edge pixels (i.e. dimension in the x-direction = readout direction) [sub-pixels]


    	double flatfieldPeak2PeakNoiseAmplitude;
    	double flatfieldWhiteNoise;
    	double flatfieldIntraPixelWidth;

    	double intraPixelSensitivity = 0.95;	// The input parameter flatfieldIntraPixelWidth is defined as the flatfield intra-pixel width at edge of pixel with 5% lower sensitivity [% of pixel size, rounded up]


    	arma::Mat<float> smearingMap;	        // Smearing map (i.e. over-scan strip)
    	arma::Mat<float> biasMap;	            // Bias map (i.e. pre-scan strip)
    	arma::Mat<float> cteMap;	            // CTE map
    	arma::Mat<float> flatfieldMap;	        // Flatfield map

    	double quantumEfficiency;	            // Quantum efficiency (in [0,1])

    	double readoutTime;                     // Readout time [s]
    	double chargeTransferTime;	            // Charge transfer time [s]

    	bool doPhotonNoise;	                    // Whether or not to apply photon noise


    	unsigned long fullWellSaturationLimit;  // Full-well saturation limit [electrons/pixel]

    	double meanCte;	                        // Mean charge-transfer efficiency
    	double readoutNoise;	                // Mean readout noise [electrons]
    	double gain;	                        // Detector gain [electrons / ADU]
    	unsigned int electronicOffset;	        // Bias or electronic offset [ADU]
    	unsigned long digitalSaturationLimit;   // Digital saturation limit [ADU / pixel]

    	string outputFilename;

    	double flatfieldSeed;
    	double readoutNoiseSeed;
    	double photonNoiseSeed;
    	double cteMapSeedRow;
    	double cteMapSeedColumn;

    	mt19937 photonNoiseGenerator;
    	mt19937 readoutNoiseGenerator;

    	poisson_distribution<int> photonNoiseDistribution;
    	normal_distribution<double> readoutNoiseDistribution;

        double internalTime;


    private:

        int imageNr;

};



#endif
