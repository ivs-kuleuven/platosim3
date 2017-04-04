#ifndef DETECTOR_H
#define DETECTOR_H

#include <string>
#include <cmath>
#include <random>
#include <functional>

#include "armadillo"

#include "Constants.h"
#include "ArrayOperations.h"
#include "Camera.h"
#include "FrontEndElectronics.h"
#include "ConfigurationParameters.h"
#include "PointSpreadFunction.h"
#include "Convolver.h"
#include "HDF5File.h"
#include "HDF5Writer.h"
#include "Logger.h"
#include "Units.h"

using namespace std;

class Camera;      // forward declaration



class Detector: public HDF5Writer 
{
    public:

    	Detector(ConfigurationParameters &configParam, HDF5File &hdf5File, Camera &camera);
    	virtual ~Detector();

    	virtual double takeExposure(int exposureNr, double startTime, double exposureTime);
    	void configure(ConfigurationParameters &configParam);

   		pair<double, double> pixelToFocalPlaneCoordinates(double row, double column);
   		pair<double, double> focalPlaneToPixelCoordinates(double xFP, double yFP);

   		pair<double, double> getFocalPlaneCoordinatesOfSubfieldCenter();
   		tuple<double, double, double, double, double, double, double, double> getFocalPlaneCoordinatesOfSubfieldCorners();

   		double getSolidAngleOfOnePixel(double plateScale);
   		double getOrientationAngle();

    	virtual tuple<bool, double, double> addFlux(double xFP, double yFP, double flux) = 0;
    	virtual void addFlux(double flux) = 0;


   		bool isInSubfield(double xFPmm, double yFPmm);


    protected:

        virtual void integrateLight(int exposureNr, double startTime, double exposureTime) = 0;

    	virtual void generateThroughputMap();
        virtual void generateGain();

    	virtual void applyFlatfield() = 0;
    	virtual void applyThroughputEfficiency();

    	virtual void readOut(float exposureTime);
    	virtual void addPhotonNoise();
    	virtual void applyFullWellSaturation();
    	virtual void applyCTI();
    	virtual void applyOpenShutterSmearing(float exposureTime);
    	virtual void addReadoutNoise();
    	virtual void applyGain();
    	virtual void addElectronicOffset();
    	virtual void applyDigitalSaturation();

    	void applySimpleCTImodel();
    	void applyShort2013CTImodel();

    	void applyParticulateContamination();
    	void applyMolecularContamination();

    	void setSubfield(const arma::Mat<float> &subfield);
    	arma::Mat<float> getSubfield();

    	virtual void initHDF5Groups() override;
    	void writePixelMapsToHDF5(int exposureNr);

        void fastForwardReadoutNoiseGeneratorToExposure(int beginExposureNr);
        void fastForwardPhotonNoiseGeneratorToExposure(int beginExposureNr);

        virtual double getTemperature();

    	arma::Mat<float> pixelMap;              // Pixel map, excl. edge pixels
    	arma::Mat<float> smearingMap;           // Smearing map (i.e. over-scan strip)
    	arma::Mat<float> biasMap;               // Bias map (i.e. pre-scan strip)
    	arma::Mat<float> throughputMap; 		// Throughput efficiency map, due to vignetting, particulate & molecular contamination, and quantum efficiency

    	unsigned int numRows; 					// Nr of rows of the detector (= size in y-direction) including non-exposed ones [pixels]
    	unsigned int numColumns; 				// Nr of columns of the detector (= size in x-direction = readout direction) [pixels]
    	unsigned int numRowsPixelMap; 			// Nr of rows in the subfield excl. edge pixels (= size the y-direction) [pixels]
    	unsigned int numColumnsPixelMap; 		// Nr of columns in the subfield excl. edge pixels (= size in the x-direction = readout direction) [pixels]
    	unsigned int numRowsSmearingMap; 		// Nr of rows in the smearing overscan strip [pixels]
    	unsigned int numRowsBiasMap; 			// Nr of rows in the bias prescan strip [pixels]

    	double originOffsetY; 					// Y-coordinate of the detector origin from the centre of the optical plane [mm]
    	double originOffsetX; 					// X-coordinate of the detector origin from the centre of the optical plane [mm]
    	unsigned int subFieldZeroPointRow; 		// Position of the subfield zeropoint w.r.t. the complete detector in the row direction [pixels]
    	unsigned int subFieldZeroPointColumn; 	// Position of the subfield zeropoint w.r.t. the complete detector in the column direction [pixels]
    	double orientationAngle; 				// Orientation angle of the detector w.r.t. the orientation of the focal plane, measured counterclockwise [radians]

    	double pixelSize;	                     // Pixel size [microns]
    	unsigned int numEdgePixels; 				 // Nr of pixels to extend the subfield on each side, to account for the edge effect

    	double polarizationEfficiency;			 // Efficiency due to polarisation at the reference angle (in [0,1])
    	double expectedValueVignetting;          // Expected value of the throughput efficiency due to vignetting (int [0,1])
    	double refAnglePolarization;			 // Reference angle for the polarisation [degrees]
    	double expectedValuePolarization;		 // Expected value of the throughput efficiency due to polarisation
    	double particulateContaminationEfficiency;	// Efficiency of particulate contamination (in [0,1])
    	double molecularContaminationEfficiency; // Efficiency of molecular contamination (in [0,1])
    	double quantumEfficiency;	             // Quantum efficiency at the reference angle (in [0,1])
    	double refAngleQuantumEfficiency;        // Reference angle for quantum efficiency [degrees]
    	double expectedValueQuantumEfficiency;   // Expected value of the throughput efficiency due to quantum efficiency
    	double readoutTime;                      // Readout time [s]
    	double readoutNoise;                     // Mean readout noise [electrons]
    double refValueGain;                     // Detector gain [µV/e-]
    double gainStability;                    // Gain stability [µV/e-]
    double gainThreeSigma;					// Allowed difference (3 sigma) in gain between the left and the right half of the detector [% of the reference value]
    double refValueGainLeft;                 // Reference value for the gain on the ACD reading the left-hand side of the detector [µV/e-]
    double refValueGainRight;                // Reference value for the gain on the ACD reading the right-hand side of the detector [µV/e-]
    	unsigned long fullWellSaturationLimit;   // Full-well saturation limit [electrons/pixel]
    	unsigned int electronicOffset;           // Bias or electronic offset [ADU]
    	unsigned long digitalSaturationLimit; 	 // Digital saturation limit [ADU / pixel]

    	string CTImodel;
    	double meanCte;               			// Mean charge-transfer efficiency  (in [0,1])
    	double beta;  							// Beta exponent in Short et al., MNRAS 430, 3078-3085 (2010).
    	double temperature;                     // Temperature of the detector
    	unsigned int numTrapSpecies; 			// Number of different trap species included in the Short2010 model
    	vector<double> trapDensity;			 	// For each trap species: the trap density [traps/pixel]
    	vector<double> trapCaptureCrossSection; // For each trap species: the trap capture cross section [m^2]
    	vector<double> releaseTime; 			// For each trap species: the electron release time [s]

    	bool includePhotonNoise;           		// Whether or not to include photon noise
    	bool includeReadoutNoise;               // Include readout noise [yes or no]
    	bool includeCTIeffects;                 // Include CTI effects [yes or no]
    	bool includeOpenShutterSmearing; 		// Include trails due reading out with an open shutter
    	bool includeQuantumEfficiency;          // Include loss of throughput due to quantum efficiency
    	bool includeVignetting;  				// Include brightness attenuation due to vignetting
    	bool includePolarization;				// Include loss of throughput due to polarisation
    	bool includeParticulateContamination;	// Include loss of throughput due to particulate contamination
    	bool includeMolecularContamination;		// Include loss of throughput due to molecular contamination
    	bool includeFullWellSaturation; 		// Whether or not full well saturation should be applied
    	bool includeDigitalSaturation; 			// Whether or not digital saturation should be applied

        int beginExposureNr;                    // Sequential number of the very first exposure. See yaml input file.

        double nominalOperatingTemperature;
    	double internalTime;

    	long readoutNoiseSeed;
    	long photonNoiseSeed;
        long gainSeed;

    	mt19937 photonNoiseGenerator;
    	mt19937 readoutNoiseGenerator;

    	poisson_distribution<long> photonNoiseDistribution;
    	normal_distribution<double> readoutNoiseDistribution;
 
        Camera &camera;
		FrontEndElectronics *frontEndElectronics;

    private:


};

#endif
