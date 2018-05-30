#include "Clock.h"

/**
 * \brief constructor for the observer class
 */
Observer::Observer()
{

}

/**
 * \brief destructor for the observer class
 */
Observer::~Observer()
{

}


/**
 * \brief constructor for the clock class
 */
Clock::Clock(string inputFilename)
{
	// clock acts as the scheduler when using the simulator parallel,
	// as does Simulation acts in the "serial" way
	// therefore it needs some input data, which it gets from the first
	// inputfile within the inputfile list

	ConfigurationParameters configParams(inputFilename);

	configure(configParams);

	endSimulation = false;
}

/**
 * \brief destructor for the clock class
 */
Clock::~Clock()
{

}

/**
 * \function that links an object to the observer list
 */
void Clock::attach(Observer* o)
{
    std::cout << "attach observer to jitter thread" << std::endl;

    observers.push_back(o);
}

/**
 * \function that removes an object to the observer list
 */
void Clock::detach(Observer* o)
{
    observers.erase(std::remove(observers.begin(), observers.end(), o), observers.end());
}

/**
 * \function that calls to update every object within the observer-list
 */
void Clock::notify()
{
    #pragma omp parallel for
    for (int i = 0; i < observers.size(); i++)
    {
        observers.at(i)->update(this, currentJitterStep);
    }
}


/**
 * \brief initial function for all simulations as a whole. when the clock is updated all 
 *  simulation objects are notified to make the next step;
 */
void Clock::startSimulation()
{
	// process simulation steps, until there 
	while(endSimulation != true)
	{
		// update the jitter
		if (!updateJitter())
		{
			// whenever there is new input, trigger the next simulation 
			// step in all observing simualtion instance
			notify();
		}	
	}
}

void Clock::updateJitter()
{
	int currentExposures = simualtionTime/(exposureTime + readoutTime);

	if (currentExposures >= numExposures)
	{
		endSimulation = true;
	}

	return endSimulation;

}

/**
 * \brief Configure the Simulation object using the input parameter file
 * 
 * \param configParams  Contains all configuration parameters from the input file
 */

void Clock::configure(ConfigurationParameters &configParams)
{
    exposureTime      = configParams.getDouble("ObservingParameters/ExposureTime"); 
    beginExposureNr   = configParams.getInteger("ObservingParameters/BeginExposureNr");
    numExposures      = configParams.getInteger("ObservingParameters/NumExposures");
    useJitter         = configParams.getBoolean("Platform/UseJitter");
    useJitterFromFile = configParams.getBoolean("Platform/UseJitterFromFile");
    includeFieldDistortion = configParams.getBoolean("Camera/IncludeFieldDistortion"); // do we want to do this or should this be asked to Camera?
    useDrift          = configParams.getBoolean("Telescope/UseDrift");  
    useDriftFromFile  = configParams.getBoolean("Telescope/UseDriftFromFile");  
    psfModel          = configParams.getString("PSF/Model");
    useFeeTemperatureFromFile = configParams.getString("FEE/Temperature") == "FromFile";
    useFeeNominalTemperature = configParams.getString("FEE/Temperature") == "Nominal";
    useDetectorTemperatureFromFile = configParams.getString("CCD/Temperature") == "FromFile";
    useDetectorNominalTemperature = configParams.getString("CCD/Temperature") == "Nominal";
    readoutTime       = configParams.getDouble("CCD/ReadoutTime"); 
}