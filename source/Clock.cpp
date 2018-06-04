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

	simulationTime = 0;

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
void Clock::attach(Observer* observerInstance)
{
    observers.push_back(observerInstance);
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
	// this uses all usable cores of the processor to speed up the simulations
    #pragma omp parallel for
    for (int i = 0; i < observers.size(); i++)
    {
        observers.at(i)->update(simulationTime);
    }

    simulationTime += timeStep;
}


/**
 * \brief initial function for all simulations as a whole. when the clock is updated all 
 *  simulation objects are notified to make the next step;
 */
void Clock::startSimulation()
{
	jitterInstance = observers.at(0)->getJitterInstance();
	driftInstance = observers.at(0)->getDriftInstance();

	// process simulation steps, until there 
	while(endSimulation != true)
	{
		if (!waitForNextStep())
		{
			// whenever there is new input, trigger the next simulation 
			// step in all observing simualtion instance

			notify();
		}	
	}

	// delete the drift- and jitter generator
    delete jitterInstance;
    delete driftInstance;
}


/**
 * \brief function that waits for new jitter values to be available 
 *        it returns false, when there are no more exposures to be processed
 */
bool Clock::waitForNextStep()
{
	timeStep = min(exposureTime, min(jitterInstance->getHeartbeatInterval(), driftInstance->getHeartbeatInterval()));

	int currentExposures = simulationTime/(exposureTime + readoutTime);

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
    numExposures      = configParams.getInteger("ObservingParameters/NumExposures");
    readoutTime       = configParams.getDouble("CCD/ReadoutTime"); 
}

