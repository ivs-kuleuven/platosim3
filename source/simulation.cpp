
#include "Simulation.h"


// Constructor

Simulation::Simulation(string inputFileName, string outputFileName)
{
    // Parse inputfile, get ConfigursationParameters objects

    hdf5File.open(outputFileName);

    // Initialise the spacecraft components

    detector   = new Detector(hdf5File);
    //camera    = new Camera(cameraConfigurationParameters);
    //telescope = new Telescope(telescopeConfigurationParameters);
    //platform  = new Platform(platformConfigurationParameters);
    //sky       = new Sky(skyConfigurationParameters);

    Nexposures = 3;        // hardcoded for the moment
    exposureTime = 22.0;   // hardcoded for the moment

}






// Destructor

Simulation::~Simulation()
{
    delete detector;
    //delete camera;
    //delete telescope;
    //delete platform;
    //delete sky;

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

