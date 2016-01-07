
#ifndef DETECTOR_H
#define DETECTOR_H

#include "timeticker.h"
#include "hdf5writer.h"
#include "camera.h"
#include "subfield.h"
#include <math.h>

using namespace std;


class Detector : public TimeTicker, Hdf5Writer
{
    public:

        Detector(ConfigurationParameters configurationParameters, Camera camera);
        virtual ~Detector();

        virtual void takeExposure(double startTime, double exposureTime);

    protected:

        virtual void reset();
    
        // Integrate light

        virtual void integrateLight(double startTime, double exposureTime);	// Integration (incl. jitter + drift) + background

        virtual void addFlux(double xCoords, double yCoords, double flux);
        virtual bool isInSubPixelMap(double row, double column);
        virtual void addFlux(double flux);

        virtual void applyFlatfield();
        virtual void rebin();

        // Read out

        virtual void readOut(double exposureTime);

        virtual void applyQuantumEfficiency(double exposureTime);
    	virtual void addPhotonNoise();
    	virtual void applyFullWellSaturation();
    	virtual void applyCte();
    	virtual void applyOpenShutterSmearing();
    	virtual void addReadoutNoise();
    	virtual void applyGain();
    	virtual void addElectronicOffset();	// Bias
    	virtual void applyDigitalSaturation();
    




        Camera camera;





        // Detector specific information

        double **pixelMap; 	// Pixel map, excl. edge pixels

    	unsigned int sizeX;// Number of columns of the detector (i.e. dimension in the x-direction = readout direction) [pixels]
    	unsigned int sizeY;// Number of rows of the detector (i.e. dimension in the y-direction) [pixels]
    	// unsigned int numSubPixelsPerPixels; // -> SubField?
    	double originOffsetX;// Offset of the detector origin from the centre of the optical plane in the x-direction [mm]
    	double originOffsetY;// Offset of the detector origin from the centre of the optical plane in the y-direction [mm]
    	double orientationAngle;// Orientation angle of the detector w.r.t. the orientation of the focal plane, measured counterclockwise [degrees]

    	unsigned int numSmearingOverscanRows;	// Number of rows in the over-scan strip
    	unsigned int numBiasPrescanRows;	// Number of rows in the pre-scan strip

    	double pixelSize;	// Pixel size [microns]
    	// double pixelScale;	// Nominal pixel scale [arcsec/pixel] // -> replace by plate scale in Camera





    	// Sub-field specific information

    	double subFieldZeroPointX;	 // Position of the sub-field zeropoint w.r.t. the complete detector in the x-direction [pixels]
    	double subFieldZeroPointY;	 // Position of the sub-field zeropoint w.r.t. the complete detector in the y-direction [pixels]

    	// Size of the sub-field in both directions [pixels]

    	int subFieldSizeX;	// Number of columns in the sub-field at pixel level and excl. edge pixels  (i.e. dimension in the x-direction = readout direction)
    	int subFieldSizeY;	// Number of rows in the sub-field at pixel leval and excl. edge pixels (i.e. dimension in the y-direction)

    	int numSubPixelsPerPixel;	// Number of sub-pixels per pixel
    	int numEdgePixels; // Number of pixels to extend the sub-field on each side, to accoutn for the edge effect

    	// Sub-pixel map and its dimensions

    	double **subPixelMap;	// Sub-pixel map, incl. edge pixels

    	int subPixelMapSizeX;	// Number of columns in the sub-field at sub-pixel level and incl. edge pixels (i.e. dimension in the x-direction = readout direction)
    	int subPixelMapSizeY;	// Number of rows in the sub-field at sub-pixel level and incl. edge pixels (i.e. dimensions in the y-direction)





    	double **smearingMap;	// Smearing map (i.e. over-scan strip)
    	double **biasRegisterMap;	// Bias register map (i.e. pre-scan strip)
    	double **cteMap;	// CTE map
    	double **flatfieldMap;	// Flatfield map

    	double quantumEfficiency;	// Quantum efficiency (in [0,1])

    	double readoutTime; // Readout time [s]
    	double chargeTransferTime;	// Charge transfer time [s]

    	bool doPhotonNoise;	// Whether or not to apply photon noise


    	unsigned long fullWellSaturationLimit;	// Full-well saturation limit [electrons/pixel]

    	//double flatfieldPeak2PeakNoise;	// Fractional peak-to-peak amplitude of the sub-pixel non-uniform sensitivity response
    	//double flatfieldSubPixelNoise;	// White-noise component of the sub-pixels
    	//double flatfieldIntraPixelWidth;	// Width of the central part of a pixel that is affected by sensitivity loss < 5% [percentage of a pixel (ceiled)]

    	double meanCte;	// Mean charge-transfer efficiency
    	double readoutNoise;	// Mean readout noise [electrons]
    	double gain;	// Detector gain [e-/ADU]
    	unsigned int electronicOffset;	// Bias or electronic offset [ADU]
    	unsigned long digitalSaturationLimit;	// Digital saturation limit [ADU/pixel]

    	string outputFilename;





    	// Seeds for random generation

    	double readoutNoiseSeed;
    	double photonNoiseSeed;
    	double flatfieldSeed;
    	// double cosmicHitSeed;
    	double cteMapSeedX;
    	double cteMapSeedY;





        double internalTime;


    private:

};



#endif
