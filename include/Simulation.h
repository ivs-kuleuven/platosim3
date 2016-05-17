#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <string>

#include "Logger.h"
#include "HDF5File.h"
#include "Detector.h"
#include "Camera.h"
#include "Telescope.h"
#include "Platform.h"
#include "Sky.h"
#include "JitterGenerator.h"
#include "JitterFromFile.h"
#include "JitterFromRedNoise.h"
#include "ConfigurationParameters.h"
#include "version.h"


using namespace std;



class Simulation
{
    public:

        Simulation(string inputFilename, string outputFilename);
        ~Simulation();
        virtual void run(double startingTime = 0.0);
        virtual void configure(ConfigurationParameters &configParams);

    protected:

        virtual void writeInputParametersToHDF5(ConfigurationParameters &configParams);
        virtual void writeVersionInformationToHDF5();

    private:

        double currentTime;
        double exposureTime;
        int Nexposures;
        bool useJitterFromFile;

        JitterGenerator *jitterGenerator;
        Platform *platform;
        Telescope *telescope;
        Sky *sky;
        Camera *camera;
        Detector *detector;

        HDF5File hdf5File;

};



#endif
