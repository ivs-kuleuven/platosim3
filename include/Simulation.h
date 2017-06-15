#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <string>

#include "Logger.h"
#include "HDF5File.h"
#include "TemperatureGenerator.h"
#include "TemperatureFromFile.h"
#include "NominalTemperature.h"
#include "Detector.h"
#include "DetectorWithMappedPSF.h"
#include "DetectorWithAnalyticGaussianPSF.h"
#include "DetectorWithAnalyticNonGaussianPSF.h"
#include "Camera.h"
#include "Telescope.h"
#include "Platform.h"
#include "Sky.h"
#include "JitterGenerator.h"
#include "JitterFromFile.h"
#include "JitterFromRedNoise.h"
#include "DriftGenerator.h"
#include "ThermoElasticDriftFromFile.h"
#include "ThermoElasticDriftFromRedNoise.h"
#include "ConfigurationParameters.h"
#include "version.h"


using namespace std;



class Simulation
{
    public:

        Simulation(string inputFilename, string outputFilename);
        ~Simulation();
        virtual void run();
        virtual void configure(ConfigurationParameters &configParams);

    protected:

        virtual void writeInputParametersToHDF5(ConfigurationParameters &configParams);
        virtual void writeVersionInformationToHDF5();
        virtual void writeStarCatalogToHDF5();

    private:

        double currentTime;
        double exposureTime;
        double readoutTime;

        int beginExposureNr;                 // sequential number of first exposure. useful for slurm parallellisation
        int numExposures;                    // Number of exposures

        bool useJitterFromFile;
        bool includeFieldDistortion;
        bool useDriftFromFile;
        bool useFeeTemperatureFromFile;
        bool useFeeNominalTemperature;
        bool useDetectorTemperatureFromFile;
        bool useDetectorNominalTemperature;
        string psfModel;

        JitterGenerator *jitterGenerator;
        DriftGenerator *driftGenerator;
        TemperatureGenerator *feeTemperatureGenerator;
        TemperatureGenerator *detectorTemperatureGenerator;
        Platform *platform;
        Telescope *telescope;
        Sky *sky;
        Camera *camera;
        Detector *detector;

        HDF5File hdf5File;

};



#endif
