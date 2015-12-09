
#include "simulation.h"


// Constructor

Simulation::Simulation(string inputFileName)
{
    // Parse inputfile, get ConfigursationParameters objects


    // Initialise the spacecraft components

    detector   = new Detector(detectorConfigurationParameters);
    camera     = new Camera(cameraConfigurationParameters);
    telescope  = new Telescope(telescopeConfigurationParameters);
    platform   = new Platform(platformConfigurationParameters);
    sky        = new Sky(skyConfigurationParameters);
    
}





// Destructor

Simulation::~Simulation()
{
    delete detector;
    delete camera;
    delete telescope;
    delete platform;
    delete sky;
}






// Simulation::run()
//
// PURPOSE:
// 
// INPUT:
//
// OUTPUT:
//

void Simulation::run(double startingTime)
{
    // Get the super-resolution subfield

    SubField subField = detector->getSubField();

    // Initialise the PSF for this particular subfield

    camera->initPsf(subField);

    // Loop over all exposures

    for (int n = 0; n < Nexposures; n++)
    {
        detector->takeExposure();
    }
}