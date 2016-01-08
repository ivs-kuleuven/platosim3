
#include "simulation.h"


// Constructor

Simulation::Simulation(string inputFileName)
{
    // Parse inputfile, get ConfigursationParameters objects


    // Initialise the spacecraft components

    detector   = new Detector();
    //camera     = new Camera(cameraConfigurationParameters);
    //telescope  = new Telescope(telescopeConfigurationParameters);
    //platform   = new Platform(platformConfigurationParameters);
    //sky        = new Sky(skyConfigurationParameters);

    Nexposures = 3;   // hardcoded for the moment
    
}





// Destructor

Simulation::~Simulation()
{
    delete detector;
    //delete camera;
    //delete telescope;
    //delete platform;
    //delete sky;
}







// Simulation::run()
//
// PURPOSE: Loop over all exposures
// 
// INPUT: 
//     . startingTime: begin time of the very first exposure. Time is expressed in seconds
//                     in the rest of the code.
//
// OUTPUT:
//     . None
//

void Simulation::run(double startingTime)
{
    currentTime = startingTime;

    // Loop over all exposures

    for (int n = 0; n < Nexposures; n++)
    {
        detector->takeExposure();
    }
}

