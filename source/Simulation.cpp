
#include "Simulation.h"


// Constructor

Simulation::Simulation(string inputFilename, string outputFilename)
{
    // Parse the configuration parameters file

    ConfigurationParameters configParams(inputFilename);

    // Open the HDF5 output file where the images will be written

    hdf5File.open(outputFilename);

    // Initialise the spacecraft components

    //platform  = new Platform(configParams);
    telescope  = new Telescope(configParams, hdf5File);
    sky        = new Sky(configParams);
    camera     = new Camera(configParams, hdf5File, *telescope, *sky);
    detector   = new Detector(configParams, hdf5File, *camera);

    Nexposures = 3;        // hardcoded for the moment
    exposureTime = 22.0;   // hardcoded for the moment

}






// Destructor

Simulation::~Simulation()
{
    // Delete order is the inverse of the order in which they were created

    delete detector;
    delete camera;
    delete telescope;
    delete sky;
    //delete platform;
    

    hdf5File.close();
}







// Simulation::run()
//
// PURPOSE: Loop over all exposures
// 
// INPUT: 
//     . startTime: begin time of the very first exposure. Time is expressed in seconds
//                  in the rest of the code.
//
// OUTPUT:
//     . None
//

void Simulation::run(double startTime)
{
    currentTime = startTime;

    // Loop over all exposures

    for (int n = 0; n < Nexposures; n++)
    {
        Log.info("Simulation: Starting exposure " + to_string(n));
        
        detector->takeExposure(startTime, exposureTime);
    }
}

