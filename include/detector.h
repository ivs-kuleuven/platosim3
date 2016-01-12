
#ifndef DETECTOR_H
#define DETECTOR_H

#include <string>
#include <cmath>

using namespace std;


class Detector
{
    public:

        Detector();
        virtual ~Detector();

        virtual void takeExposure(double startTime, double exposureTime);

    protected:

        virtual void reset();
    
        // Integrate light

        virtual void integrateLight(double startTime, double exposureTime);

        virtual void addFlux(double xCoords, double yCoords, double flux);
        virtual bool isInSubPixelMap(double row, double column);
//        virtual void addFlux(double flux);

        virtual void applyFlatfield();
        virtual void rebin();

        // Read out

        virtual void readOut(double exposureTime);

        virtual void applyQuantumEfficiency();
        virtual void addSkyBackground(double exposureTime);
    	virtual void addPhotonNoise();
    	virtual void applyFullWellSaturation();
    	virtual void applyCte();
    	virtual void applyOpenShutterSmearing();
    	virtual void addReadoutNoise();
    	virtual void applyGain();
    	virtual void addElectronicOffset();	// Bias
    	virtual void applyDigitalSaturation();
    


        // Detector specific information

        double **pixelMap; 	                   // Pixel map, excl. edge pixels

        unsigned int numRows;                  // Number of rows of the detector (i.e. dimension in the y-direction) [pixels]
    	unsigned int numColumns;               // Number of columns of the detector (i.e. dimension in the x-direction = readout direction) [pixels]
    	double originOffsetRow;                // Offset of the detector origin from the centre of the optical plane in the row direction [mm]
    	double originOffsetColumn;             // Offset of the detector origin from the centre of the optical plane in the column direction [mm]
    	double orientationAngle;               // Orientation angle of the detector w.r.t. the orientation of the focal plane, measured counterclockwise [degrees]
 
    	unsigned int numRowsSmearingMap;       // Number of rows in the smearing over-scan strip [pixels]
    	unsigned int numRowsBiasMap;	       // Number of rows in the bias pre-scan strip [pixels]

    	double pixelSize;	                   // Pixel size [microns]
    

    	// Sub-field specific information

    	double subFieldZeroPointRow;	       // Position of the sub-field zeropoint w.r.t. the complete detector in the row direction [pixels]
    	double subFieldZeroPointColumn;	       // Position of the sub-field zeropoint w.r.t. the complete detector in the column direction [pixels]
    	int numRowsSubField;	               // Number of rows in the sub-field at pixel level and excl. edge pixels (i.e. dimension in the y-direction) [pixels]
    	int numColumnsSubField;	               // Number of columns in the sub-field at pixel level and excl. edge pixels  (i.e. dimension in the x-direction = readout direction) [pixels]
    	int numSubPixelsPerPixel;	           // Number of sub-pixels per pixel
    	int numEdgePixels;                     // Number of pixels to extend the sub-field on each side, to account for the edge effect


    	// Sub-pixel map and its dimensions

    	double **subPixelMap;	               // Sub-pixel map, incl. edge pixels

    	int numRowsSubPixelMap;	               // Number of rows in the sub-field at sub-pixel level and incl. edge pixels (i.e. dimensions in the y-direction) [sub-pixels]
    	int numColumnsSubPixelMap;	           // Number of columns in the sub-field at sub-pixel level and incl. edge pixels (i.e. dimension in the x-direction = readout direction) [sub-pixels]




    	double **smearingMap;	               // Smearing map (i.e. over-scan strip)
    	double **biasMap;	                   // Bias map (i.e. pre-scan strip)
    	double **cteMap;	                   // CTE map
    	double **flatfieldMap;	               // Flatfield map

    	double quantumEfficiency;	           // Quantum efficiency (in [0,1])

    	double readoutTime;                    // Readout time [s]
    	double chargeTransferTime;	           // Charge transfer time [s]

    	bool doPhotonNoise;	                   // Whether or not to apply photon noise


    	unsigned long fullWellSaturationLimit; // Full-well saturation limit [electrons/pixel]

    	double meanCte;	                       // Mean charge-transfer efficiency
    	double readoutNoise;	               // Mean readout noise [electrons]
    	double gain;	                       // Detector gain [electrons / ADU]
    	unsigned int electronicOffset;	       // Bias or electronic offset [ADU]
    	unsigned long digitalSaturationLimit;  // Digital saturation limit [ADU / pixel]

    	string outputFilename;


    	// Seeds for random generation

    	double readoutNoiseSeed;
    	double photonNoiseSeed;
    	double flatfieldSeed;
    	double cteMapSeedRow;
    	double cteMapSeedColumn;





        double internalTime;


    private:

};



#endif
