#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <string>

#include "Logger.h"
#include "HDF5File.h"
#include "Detector.h"
#include "Camera.h"
#include "Telescope.h"
//#include "Platform.h"
#include "Sky.h"
#include "ConfigurationParameters.h"


using namespace std;



class Simulation
{
    public:

        Simulation(string inputFilename, string outputFilename);
        ~Simulation();
        virtual void run(double startingTime = 0.0);

    protected:

        
    private:

        double currentTime;
        double exposureTime;
        int Nexposures;


        Detector  *detector;
        Camera    *camera;
        Telescope *telescope;
        //Platform  *platform;
        Sky       *sky;

        HDF5File hdf5File;

};



#endif
