
#include "controller.h"


// Constructor

Controller::Controller(string inputFileName)
{
    // Parse inputfile, get ConfigruationParameters objects


    // Initialise the spacecraft components

    detector   = new Detector(detectorConfigurationParameters);
    camera     = new Camera(cameraConfigurationParameters);
    telescope  = new Telescope(telescopeConfigurationParameters);
    platform   = new Platform(platformConfigurationParameters);
    sky        = new Sky(skyConfigurationParameters);
    
}





// Destructor

Controller::~Controller()
{
    delete detector;
    delete camera;
    delete telescope;
    delete platform;
    delete sky;
}






// Controller::run()
//
// PURPOSE:
// 
// INPUT:
//
// OUTPUT:
//

void Controller::run(double startingTime)
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