#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <string>
#include "logger.h"
#include "hdf5file.h"
#include "detector.h"
//#include "camera.h"
//#include "telescope.h"
//#include "platform.h"
//#include "sky.h"
//#include "configurationparameters.h"


using namespace std;



class Simulation
{
    public:

        Simulation(string inputFileName, string outputFileName);
        ~Simulation();
        virtual void run(double startingTime = 0.0);

    protected:

        
    private:

        double currentTime;
        double exposureTime;
        int Nexposures;


        Detector  *detector;
        //Camera    *camera;
        //Telescope *telescope;
        //Platform  *platform;
        //Sky       *sky;

        HDF5File hdf5File;

};



#endif
