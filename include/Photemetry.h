#ifndef PHOTOMETRY_H
#define PHOTOMETRY_H

#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include <functional>
#include <complex>
#include <string>

#include "armadillo"

#include "Constants.h"
#include "Units.h"
#include "Detector.h"
#include "Camera.h"
#include "Parameter.h"

using namespace std;



class DetectorWithAnalyticNonGaussianPSF: public Detector 
{
 public:

  // Mapped Zemax PSF
  
  DetectorWithMappedPSF(ConfigurationParameters &configParam,
			HDF5File &hdf5file,
			Camera &camera,
			TemperatureGenerator &feeTemperatureGenerator,
			TemperatureGenerator &detectorTemperatureGenerator,
			double readoutTimeBeforeNextExposure,
			double readoutTimeDuringNextExposure);
  ~DetectorWithMappedPSF();

  // Analytic PSF
  
  DetectorWithAnalyticNonGaussianPSF(ConfigurationParameters &configParam,
				     HDF5File &hdf5File,
				     Camera &camera,
				     TemperatureGenerator &feeTemperatureGenerator,
				     TemperatureGenerator &detectorTemperatureGenerator,
				     double readoutTimeBeforeNextExposure,
				     double readoutTimeDuringNextExposure);
  ~DetectorWithAnalyticNonGaussianPSF();

 protected:

	// Photometry
  
        bool includePhotometry;             // Whether or not to include on-the-fly photometry
        int contaminationRadius;            // Stars outside the radius are never considered a contaminant of the main target [pix] 

	// Aperture photometry

	double maskUpdateInterval;          // All photometric masks will be updated every xx days  [days]
        vector<unsigned int> photStarIDs;   // Star IDs for which you need photometry

        // The following maps contain for each star ID, a vector with the photometry information for each exposure.
        // E.g. maskSizeTarget[1234][10] contains the nr of pixels of the mask of target 1234 for exposure 10.

        map<unsigned int, vector<unsigned int>> exposureNrOfMaskUpdate;   // Photometric mask is updated once in a while
        map<unsigned int, vector<unsigned int>> maskSizeTarget;           // Nr of pixels within the mask for each target
        map<unsigned int, vector<double>> inputFluxTarget;                // Flux of the target as computed from the (variable) Vmag
        map<unsigned int, vector<double>> estimatedFluxTarget;            // Estimated flux for each target
        map<unsigned int, vector<double>> varFluxTarget;                  //  Variance of the flux for each target
        map<unsigned int, vector<double>> NSRtarget;                      // Noise/Signal ratio of the flux of each target
    
        // The indices (row and column) of the masks are stored as [starID][exposureNr][index]
  
        map<unsigned int, map<unsigned int, vector<unsigned int>>> rowIndexOfMaskOfTarget;  
        map<unsigned int, map<unsigned int, vector<unsigned int>>> colIndexOfMaskOfTarget;

    private:
  
};

#endif
