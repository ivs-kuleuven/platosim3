#ifndef CAMERA_H
#define CAMERA_H

#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <array>
#include <tuple>

#include "armadillo"

#include "ArrayOperations.h"
#include "ConfigurationParameters.h"
#include "Constants.h"
#include "Detector.h"
#include "HDF5File.h"
#include "HDF5Writer.h"
#include "Heartbeat.h"
#include "Logger.h"
#include "Exceptions.h"
#include "Polynomial1D.h"
#include "Sky.h"
#include "StringUtilities.h"
#include "Platform.h"
#include "Telescope.h"
#include "Units.h"
#include "Parameter.h"


using namespace std;


class Detector;  // forward declaration

typedef map<unsigned int, array<double, 6+1>>::iterator starInfoIterator;


class Camera : public HDF5Writer
{
    public:

        Camera(ConfigurationParameters &configParam, HDF5File &hdf5File, Platform &platform, Telescope &telescope, Sky &sky);
        virtual ~Camera();

        virtual void configure(ConfigurationParameters &configParam);
        virtual pair<double, double> exposeDetector(Detector &detector, double startTime, double exposureTime, double readoutTimeBeforeNextExposure, int binnumber, int subsubfieldx, int subsubfieldy);	//%% changed for spectral dependence, added the wavelength bin, also added subsubfield for splitting, made pair of double
	virtual void SkyBackground(Detector &detector, double startTime, double exposureTime, double readoutTimeBeforeNextExposure, int binnumber, double centerRA, double centerDec, int subsubfieldx, int subsubfieldy);  //%% added for spectral dependence, apply sky BG separately
        virtual void updateParameters(double time);

        virtual void initHDF5Groups() override;
        virtual void flushOutput() override;

        pair<double, double> skyToFocalPlaneCoordinates(double raStar, double decStar, bool useInitialOrientation=false);
        pair<double, double> focalPlaneToSkyCoordinates(double xFP, double yFP, bool useInitialOrientation=false);

        pair<double, double> undistortedToDistortedFocalPlaneCoordinates(double xFPmm, double yFPmm);
        pair<double, double> distortedToUndistortedFocalPlaneCoordinates(double xFPdist, double yFPdist);

        double getGnomonicRadialDistanceFromOpticalAxis(double xFP, double yFP);

        set<unsigned int> getAllStarIDs();

        tuple<double, double, double, double, double, double> getInfoForTheMostRecentExposureForStar(int starID);
        pair<starInfoIterator, starInfoIterator> getInfoForTheMostRecentExposureForAllStars();

        vector<double> getTotalSkyBackground();  //%% changed for spectral dependence, is now a vector with n wavelength bins

        double getFocalLength();


    protected:

        int beginExposureNr;                 // sequential number of first exposure. useful for slurm parallellisation
        int numExposures;                    // Number of exposures

        Platform &platform;
        Telescope &telescope;
        Sky &sky;

        Parameter<double> *focalLength;       // [mm]
        Parameter<double> *focalPlaneAngle;   // Orientation of the focal plane, as an angle around the optical axis  [rad] 
        Parameter<double, 3> *distortionCoef; // distortion coefficients to map undistorted to distorted coordinates.
        Parameter<double, 3> *inverseDistortionCoef; // inverse distortion coefficient to map distorted to undistorted coordinates.

        string distortionModel;               // The model used to compute the distortion  
        double plateScale;                    // [arcsec/micron]
        double throughputLambdaC;             // Central wavelength of the throughput passband [nm]

	int wavelengthBins;  //%% Added for spectral dependence, number of wavelength bins
	double binWidth;  //%% width of one wavelength bin [nm]
        double binOrigin;  //%% lower edge of first wavelength bin [nm]
	int numsubsubfieldsx;  //%% Number of subsubfields to execute separately, x
	int numsubsubfieldsy;  //%% Number of subsubfields to execute separately, x

        double internalTime;

        bool includeAberrationCorrection; // Whether or not (differential) aberration correction should be included
        string aberrationCorrectionType;  // [differential or absolute]

        bool includeFieldDistortion;      // Whether or not field distortion should be included

        double userGivenSkyBackground;    // User-set zodiacal + stellar sky background.                          [phot/pix/s]
                                          // If negative, computed by the Sky class
        double fluxOfV0Star;              // Photon flux of a V=0 (G2V) star                                      [phot/s/m^2/nm]

        double raSun;                     // Right ascension of the direction of the sun shield during the run    [rad]
        double decSun;                    // Declination of the direction of the sun shield during the run        [rad]

        bool writeStarPositions;          // Whether or not the star positions should be written to the output HDF5 file

        // detectedStarInfo[startTime][starID] contains the values (xFPmean, yFPmean, rowPixMean, colPixmean, sumFlux, Ndetections) + temp

        map<double, map<unsigned int, array<double, 6+1>>> detectedStarInfo;	//%% +1 for spetral dependence, as temperature is also tranferred
        vector<double> skyBackgroundValues;
        vector<double> transmissionEfficiencyValues;
        double totalSkyBackground;          // Total sky background [photons / pixel / exposure]

    private:

        
};



#endif
