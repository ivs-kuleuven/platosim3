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
void Clock::startSimulation(std::condition_variable* cond_var, bool* notified, bool* newStep, std::mutex* m)
{
	jitterInstance = observers.at(0)->getJitterInstance();
	driftInstance = observers.at(0)->getDriftInstance();

	// process simulation steps, until there 
	while(!endSimulation)
	{
		
		// check whether jitter is coming from a tcp connection
		if (jitterInstance->isClient())
		{
			*notified = true;
			
			// notify the tcp connection thread
			cond_var->notify_one();		

			// declare a lock for this thread
			std::unique_lock<std::mutex> lock(*m);

	        // wait for the tcp connection thread to notify this thread
	        while(!*newStep)
	        {    	
	            cond_var->wait(lock);
	        }

	        // reset to the initial parameters
	        *newStep = false;
			lock.unlock();
		}
		
		// get the next time step for the simulation
		timeStep = min(exposureTime, min(jitterInstance->getHeartbeatInterval(), driftInstance->getHeartbeatInterval()));

		// check, whether the intended number of imagettes is already reached
		int currentExposures = simulationTime/(exposureTime + readoutTime);

		if (currentExposures >= numExposures)
		{
			endSimulation = true;
		}

		// if there are still jitter steps coming from the server, keep the simulation going
		if (jitterInstance->isClient())
		{
			if(!jitterInstance->simulationEnd())
			{
				endSimulation = false;
			}
			else
			{
				endSimulation = true;	
			}
		}

		if (!endSimulation)
		{
			// whenever there is new input, trigger the next simulation 
			// step in all observing simualtion instances

			notify();
		}	
	}

	// delete the drift- and jitter generator
    delete jitterInstance;
    delete driftInstance;
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