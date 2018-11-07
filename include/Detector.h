#ifndef DETECTOR_H
#define DETECTOR_H

#include <string>
#include <cmath>
#include <random>
#include <functional>
#include <valarray>

#include "armadillo"

#include "Faddeeva.hh"

#include "Constants.h"
#include "ArrayOperations.h"
#include "Mathematics.h"
#include "Camera.h"
#include "FrontEndElectronics.h"
#include "TemperatureGenerator.h"
#include "ConfigurationParameters.h"
#include "PointSpreadFunction.h"
#include "Convolver.h"
#include "HDF5File.h"
#include "HDF5Writer.h"
#include "Logger.h"
#include "Units.h"

using namespace std;

class Camera;      // forward declaration



class IntegralOfAnalyticSignalResponse
{
    public:

        IntegralOfAnalyticSignalResponse() : size(0), n(0.) {};
        IntegralOfAnalyticSignalResponse(size_t s, double d = 0.) : size(s), n(0.), dsigma(d) {}
        virtual ~IntegralOfAnalyticSignalResponse(){};
        IntegralOfAnalyticSignalResponse& addPart(double, double, double, double, double = 0., double = 0., double = 0.);
        double operator()(unsigned, unsigned, bool = true);

    private:

        size_t size;                              // number of (sub)pixels in one dimension
        double n;                                 // normalization factor
        double dsigma;                            // Gaussian diffusion kernel width
        vector<valarray<double>> erfxr;           // evaluated error functions for x
        vector<valarray<double>> erfyr;           // evaluated error functions for y
        vector<valarray<complex<double>>> erfxc;  // evaluated complex error functions for x
        vector<valarray<complex<double>>> erfyc;  // evaluated complex error functions for y
};










class Detector: public HDF5Writer 
{
    public:

        Detector(ConfigurationParameters &configParam, HDF5File &hdf5File, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure);
        virtual ~Detector();

        virtual double takeExposure(int exposureNr, double startTime, double exposureTime);
        virtual void updateParameters(double time);
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

        double getReadoutTimeBeforeNextExposure();


    protected:

        virtual void integrateLight(int exposureNr, double startTime, double exposureTime) = 0;

        virtual void generateThroughputMap();
        virtual void checkGain();
        virtual void generateGuyonnetCoefficients();

        virtual void applyFlatfield() = 0;
        virtual void applyThroughputEfficiency();
        virtual void applyBFE();
        virtual void addDarkSignal(float exposureTime);

        virtual void readOut(float exposureTime);
        virtual void addPhotonNoise();
        virtual void addCosmics(float exposureTime);
        virtual void addCosmics(float exposureTime, arma::Mat<float> &map, int numRows, int numColumns, string area);
        virtual void applyFullWellSaturation();
        virtual void applyCTI();
        virtual void applyOpenShutterSmearing(float exposureTime);
        virtual void addReadoutNoise();
        virtual void applyQuantisation();
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

        virtual double getTemperature();

        arma::Mat<float> pixelMap;               // Pixel map, excl. edge pixels
        arma::Mat<float> smearingMap;            // Smearing map (i.e. over-scan strip)
        arma::Mat<float> biasMapLeft;            // Bias map (i.e. pre-scan strip) for the left detector half
        arma::Mat<float> biasMapRight;            // Bias map (i.e. pre-scan strip) for the right detector half
        arma::Mat<float> throughputMap;          // Throughput efficiency map, due to vignetting, particulate & molecular contamination, and quantum efficiency

        unsigned int numRows;                    // Nr of rows of the detector (= size in y-direction) including non-exposed ones [pixels]
        unsigned int numColumns;                 // Nr of columns of the detector (= size in x-direction = readout direction) [pixels]
        unsigned int numRowsPixelMap;            // Nr of rows in the subfield excl. edge pixels (= size the y-direction) [pixels]
        unsigned int numColumnsPixelMap;         // Nr of columns in the subfield excl. edge pixels (= size in the x-direction = readout direction) [pixels]
        unsigned int numRowsSmearingMap;         // Nr of rows in the smearing over-scan strip [pixels]
        unsigned int numRowsBiasMap;             // Nr of rows in the bias pre-scan strip [pixels]
        unsigned int numColumnsBiasMap;          // Nr of columns in the bias pre-scan strip [pixels]

        unsigned int firstRowExposed;            // Index of the first row that is exposed to light (different for the Fast and Normal camera's) [pixels]

        double originOffsetY;                    // Y-coordinate of the detector origin from the centre of the optical plane [mm]
        double originOffsetX;                    // X-coordinate of the detector origin from the centre of the optical plane [mm]
        unsigned int subFieldZeroPointRow;       // Position of the subfield zeropoint w.r.t. the complete detector in the row direction [pixels]
        unsigned int subFieldZeroPointColumn;    // Position of the subfield zeropoint w.r.t. the complete detector in the column direction [pixels]
        double orientationAngle;                 // Orientation angle of the detector w.r.t. the orientation of the focal plane, measured counterclockwise [radians]

        double pixelSize;                        // Pixel size [microns]
        unsigned int numEdgePixels;              // Nr of pixels to extend the subfield on each side, to account for the edge effect

        arma::Cube<float> guyonnetCoefficients;  // Coefficients a^X_ij for the BFE in Sect. 6.1 in Guyonnet et al. 2015
        double p0BFE;        					 // Value for p0 parameter in Eq. (18) in Guyonnet et al. 2015
        double p1BFE;						     // Value for p1 parameter in Eq. (18) in Guyonnet et al. 2015
        int rangeBFE;							 // How far pixels can be apart and still influence each other [pixels] (use window with dimensions 2 * range + 1)
        double refFluxBFE;                       // Reference flux for the p0 and p1 parameters for BFE [e-]


        bool includeCosmicsInSubField;           // Whether or not to include cosmic hits in the subfield
        bool includeCosmicsInSmearingMap;        // Whether or not to include cosmic hits in the (physical) overscan region
        bool includeCosmicsInBiasMap;            // Whether or not to include cosmic hits in the (virtual) prescan region
        double cosmicHitRate;					 // Cosmic hit rate [events / cm^2 / s]
        vector<double> cosmicTrailLength;		 // Interval of the length of the cosmic trails [pixels]
        vector<double> cosmicIntensity; 		 // Interval of the intensity of the cosmic trails [e-]
        double expectedValueVignetting;          // Expected value of the throughput efficiency due to vignetting (int [0,1])
        double expectedValuePolarization;        // Expected value of the throughput efficiency due to polarisation
        double particulateContaminationEfficiency;  // Efficiency of particulate contamination (in [0,1])
        double molecularContaminationEfficiency;    // Efficiency of molecular contamination (in [0,1])
        double meanQE;							 // Mean QE (over all wavelengths)
        double meanAngleDependencyQE;			 // Mean (over all pixels) of the relative efficiency due to the angle dependency of the QE
        double serialTransferTime;				 // Time to shift the content of the readout register by one pixel [s]
        double parallelTransferTime;			 // Time to shift the charges one row down in case the readout register will be read out [s]
        double parallelTransferTimeFast;	     // Time to shift the charges one row down in case the readout register will not be read out [s]
        bool isFastCamera;                       // Indicates whether or not the camera is a fast camera
        int firstRowPartialReadout;			     // First row that will be read out by the FEE in partial readout mode
        int numRowsPartialReadout;			     // Number of rows that will be read out by the FEE, starting at firstRowReadout, in partial readout mode
        double readoutNoise;                     // Mean readout noise [electrons]
        double refValueGainLeft;                 // Reference value for the gain on the ACD reading the left-hand side of the detector [µV/e-]
        double refValueGainRight;                // Reference value for the gain on the ACD reading the right-hand side of the detector [µV/e-]
        double gainStability;                    // Gain stability [µV/e-]
        double gainAllowedDifference;            // Allowed difference in gain between the left and the right half of the detector [% of the reference values]
        unsigned long fullWellSaturationLimit;   // Full-well saturation limit [electrons/pixel]
        unsigned int electronicOffset;           // Bias or electronic offset [ADU]
        unsigned long digitalSaturationLimit;    // Digital saturation limit [ADU / pixel]
        double darkCurrent;						 // Dark current [e- / s]
        double dsnu;							 // Dark signal non-uniformity
        double darkCurrentStability;             // Temperature stability of the dark current [e / K / s]

        string CTImodel;
        double meanCte;                          // Mean charge-transfer efficiency  (in [0,1])
        double beta;                             // Beta exponent in Short et al., MNRAS 430, 3078-3085 (2010).
        double temperature;                      // Temperature of the detector
        unsigned int numTrapSpecies;             // Number of different trap species included in the Short2010 model
        vector<double> trapDensity;              // For each trap species: the trap density [traps/pixel]
        vector<double> trapCaptureCrossSection;  // For each trap species: the trap capture cross section [m^2]
        vector<double> releaseTime;              // For each trap species: the electron release time [s]

        string readoutMode;                      // Readout mode (Nominal / Partial)
        double readoutTimeBeforeNextExposure;    // Duration of the readout before the next exposure can start [s]
        double readoutTimeDuringNextExposure;    // Duration of the readout when the next exposure has already started [s]

        bool includeBFE;						 // Whether or not to include the BFE
        bool includeDarkSignal;	      			 // Whether or not to include dark
        bool includePhotonNoise;                 // Whether or not to include photon noise
        bool includeReadoutNoise;                // Include readout noise [yes or no]
        bool includeCTIeffects;                  // Include CTI effects [yes or no]
        bool includeOpenShutterSmearing;         // Include trails due reading out with an open shutter
        bool includeQuantumEfficiency;           // Include loss of throughput due to quantum efficiency
        bool includeVignetting;                  // Include brightness attenuation due to vignetting
        bool includePolarization;                // Include loss of throughput due to polarisation
        bool includeParticulateContamination;    // Include loss of throughput due to particulate contamination
        bool includeMolecularContamination;      // Include loss of throughput due to molecular contamination
        bool includeFullWellSaturation;          // Whether or not full well saturation should be applied
        bool includeDigitalSaturation;           // Whether or not digital saturation should be applied
        bool includeQuantisation;                // Whether or not to include quantisation

        int beginExposureNr;                     // Sequential number of the very first exposure. See yaml input file.

        double nominalOperatingTemperature;
        double internalTime;

        long darkSignalSeed;
        long readoutNoiseSeed;
        long photonNoiseSeed;
        long cosmicSeed;

        mt19937 darkSignalGenerator;
        mt19937 darkNoiseGenerator;
        mt19937 photonNoiseGenerator;
        mt19937 readoutNoiseGenerator;
        mt19937 cosmicHitRateGenerator;
        mt19937 cosmicEntryRowGenerator;
        mt19937 cosmicEntryColumnGenerator;
        mt19937 cosmicEntryAngleGenerator;
        mt19937 cosmicTrailLengthGenerator;
        mt19937 cosmicIntensityGenerator;
        mt19937 decimalNumCosmicHitsGenerator;

        normal_distribution<double> darkSignalDistribution;
        normal_distribution<double> darkNoiseDistribution;
        poisson_distribution<long> photonNoiseDistribution;
        normal_distribution<double> readoutNoiseDistribution;
        poisson_distribution<long> cosmicHitRateDistribution;
        uniform_real_distribution<double> cosmicEntryRowDistribution;
        uniform_real_distribution<double> cosmicEntryColumnDistribution;
        uniform_real_distribution<double> cosmicEntryAngleDistribution;
        uniform_real_distribution<double> cosmicTrailLengthDistribution;
        uniform_real_distribution<double> cosmicIntensityDistribution;
        uniform_real_distribution<double> decimalNumCosmicHitsDistribution;
 
        Camera &camera;
        FrontEndElectronics *frontEndElectronics;

    private:

        TemperatureGenerator &temperatureGenerator;

};

#endif
