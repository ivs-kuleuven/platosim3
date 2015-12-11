
#ifndef DETECTOR_H
#define DETECTOR_H

#include <string>
#include "telescope.h"
#include "subfield.h"


using namespace std;


class Detector
{
    public:

        Detector(ConfigurationParameters configurationParameters, Camera camera);
        virtual ~Detector();

        virtual void takeExposure(double startTime, double exposureTime);
        SubField getSubField();

    protected:
    
        // Integrate light

        virtual void integrateLight(double startTime, double exposureTime);	// Integration (incl. jitter + drift) + background

        virtual void applyFlatfield();

        // Read out

        virtual void readOut();


        virtual void applyQuantumEfficiency();
    	virtual void addPhotonNoise();
    	virtual void applyFullWellSaturation();
    	virtual void applyCte();
    	virtual void applyOpenShutterSmearing();
    	virtual void addReadoutNoise();
    	virtual void applyGain();
    	virtual void addElectronicOffset();	// Bias
    	virtual void applyDigitalSaturation();
    




        Camera camera;
        SubField subField;

    	unsigned int sizeX;// Number of columns of the detector (i.e. dimension in the x-direction = readout direction) [pixels]
    	unsigned int sizeY;// Number of rows of the detector (i.e. dimension in the y-direction) [pixels]
    	// unsigned int numSubPixelsPerPixels; // -> SubField?
    	double originOffsetX;// Offset of the detector zeropoint from the centre of the optical plane in the x-direction [mm]
    	double originOffsetY;// Offset of the detector zeropoint from the centre of the optical plane in the y-direction [mm]
    	double orientationAngle;// Orientation angle of the detector w.r.t. the orientation of the focal plane, measured counterclockwise [degrees]

    	unsigned int numSmearingOverscanRows;	// Number of rows in the over-scan strip
    	unsigned int numBiasPrescanRows;	// Number of rows in the pre-scan strip

    	unsigned double pixelSize;	// Pixel size [microns]
    	// double pixelScale;	// Nominal pixel scale [arcsec/pixel] // -> replace by plate scale in Camera


    	double smearingMap[][];	// Smearing map (i.e. over-scan strip)
    	double biasRegisterMap[][];	// Bias register map (i.e. pre-scan strip)
    	double cteMap[][];	// CTE map
    	double flatfieldMap;	// Flatfield map

    	unsigned double quantumEfficiency;	// Quantum efficiency (in [0,1])

    	unsigned double readoutTime; // Readout time [s]
    	unsigned double chargeTransferTime;	// Charge transfer time [s]

    	bool doPhotonNoise;	// Whether or not to apply photon noise


    	unsigned long fullWellSaturationLimit;	// Full-well saturation limit [electrons/pixel]

    	//double flatfieldPeak2PeakNoise;	// Fractional peak-to-peak amplitude of the sub-pixel non-uniform sensitivity response
    	//double flatfieldSubPixelNoise;	// White-noise component of the sub-pixels
    	//double flatfieldIntraPixelWidth;	// Width of the central part of a pixel that is affected by sensitivity loss < 5% [percentage of a pixel (ceiled)]

    	unsigned double meanCte;	// Mean charge-transfer efficiency
    	unsigned double readoutNoise;	// Mean readout noise [electrons]
    	unsigned double gain;	// Detector gain [e-/ADU]
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



        double **pixelMap;                       // Par of a CCD image. Nrows x Ncols
 
        double internalTime;


    private:

};



#endif
