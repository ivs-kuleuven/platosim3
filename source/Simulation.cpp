
#include "Simulation.h"


// Constructor

Simulation::Simulation(string inputFilename, string outputFilename)
{
    // Parse the configuration parameters file

    Log.info("Simulation: reading the input parameters file");

    ConfigurationParameters configParams(inputFilename);

    // Check if the output HDF5 filename already exists. If so, complain.

    if (fileExists(outputFilename))
    {
        Log.error("Simulation: Output file name already exists. Aborting.");
        exit(1);
    }

    // Open the HDF5 output file where the images will be written

    hdf5File.open(outputFilename);

    // Configure the Simulation object using the configuration parameters file

    configure(configParams);

    // Depending on what the user requested, define the proper jitter generator

    if (useJitterFromFile)
    {
        jitterGenerator = new JitterFromFile(configParams);
    }
    else
    {
        jitterGenerator = new JitterFromRedNoise(configParams);
    }

    // Initialise the spacecraft components

    platform   = new Platform(configParams, hdf5File, *jitterGenerator);
    telescope  = new Telescope(configParams, hdf5File, *platform);
    sky        = new Sky(configParams);
    camera     = new Camera(configParams, hdf5File, *telescope, *sky);
    detector   = new Detector(configParams, hdf5File, *camera);

}






// Destructor

Simulation::~Simulation()
{
    // Delete order is the inverse of the order in which they were created

    delete detector;
    delete camera;
    delete telescope;
    delete sky;
    delete platform;
    delete jitterGenerator;
    

    hdf5File.close();
}








/**
 * \brief Configure the Simulation object using the input parameter file
 * 
 * \param configParams  Contains all configuration parameters from the input file
 */

void Simulation::configure(ConfigurationParameters &configParams)
{
    exposureTime      = configParams.getDouble("ObservingParameters/ExposureTime"); 
    Nexposures        = configParams.getInteger("ObservingParameters/NumExposures"); 
    useJitterFromFile = configParams.getInteger("Platform/UseJitterFromFile");
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

