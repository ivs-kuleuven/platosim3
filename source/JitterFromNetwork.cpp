#include "JitterFromNetwork.h"



/**
 * \brief Constructor of the jitter from network class
 * 
 *       
 */
JitterFromNetwork::JitterFromNetwork(ConfigurationParameters &configParams, double readoutTimeBeforeNextExposure, std::condition_variable* cond_var, std::mutex* m, bool* notified, bool* newStep)
{
    // Set the configuration parameters

    Log.info("JitterFromNetwork: object created");
  
    configure(configParams);

    condVarPointer = cond_var;
    mutexPointer = m;
    notifiedPointer = notified;
    newStepPointer = newStep;

    heartbeatInterval = 0.01;

}


/**
 * \brief Destructor
 */

JitterFromNetwork::~JitterFromNetwork()
{

}


/**
 * \brief Configure this object using the parameters from the input parameters file
 * 
 * \param configParams  The configuration parameters
 */

void JitterFromNetwork::configure(ConfigurationParameters &configParams)
{
	oldTimeStep = -0.1;
    oldYaw = 0.0;
    oldPitch = 0.0;
    oldRoll = 0.0;

    currentTimeStep = 0.0;
    currentYaw = 0.0;
    currentPitch = 0.0;
    currentRoll = 0.0;
}



/**
 * \brief returns the jitter step (yaw, pitch, roll) at the given time value
 * 
 * \
 */
tuple<double, double, double> JitterFromNetwork::getNextYawPitchRoll(double time)
{
	// ask for the server thread to get a new jitter step until the new step has a larger time stamp than time

	while(time >= currentTimeStep)
	{
		Log.info("JitterFromNetwork: time: " + to_string(time));
		Log.info("JitterFromNetwork: currentTime: " + to_string(currentTimeStep));

		*notifiedPointer = true;

		*newStepPointer = false;
			
		Log.info("JitterFromNetwork: notify jitter thread");

		// notify the tcp connection thread
		condVarPointer->notify_one();	

		// declare a lock for this thread
		std::unique_lock<std::mutex> lock(*mutexPointer);

		Log.info("JitterFromNetwork: wait for new step");

	    // wait for the tcp connection thread to notify this thread
	    while(!*newStepPointer)
	    {    	
	        condVarPointer->wait(lock);
	    }

		Log.info("JitterFromNetwork: notification from jitter thread received");

	    if (abs(time - currentTimeStep) < 0.000001)
	    {
	        Log.info("JitterFromNetwork: time: " + to_string(currentTimeStep) + "; yaw: " + to_string(currentYaw*1000) + "x10^-3; pitch: " + to_string(currentPitch*1000) + "x10^-3; roll: " + to_string(currentRoll*1000) + "x10^-3");

			return make_tuple(currentYaw, currentPitch, currentRoll);
	    }
	    else if (currentTimeStep > time)
	    {
	        // interpolate the jitterstep dependent on the old and new jitterstep and the currentTime

	        const double weight1 = (time - oldTimeStep) / (currentTimeStep - oldTimeStep);
	        const double weight2 = (currentTimeStep - time) / (currentTimeStep - oldTimeStep);
	        const double newYaw   = oldYaw   * weight2 + currentYaw   * weight1;
	        const double newPitch = oldPitch * weight2 + currentPitch * weight1;
	        const double newRoll  = oldRoll  * weight2 + currentRoll  * weight1;

	        Log.info("JitterFromNetwork: time: " + to_string(currentTimeStep) + "; yaw: " + to_string(newYaw*1000) + "x10^-3; pitch: " + to_string(newPitch*1000)+ "x10^-3; roll: " + to_string(newRoll*1000) + "x10^-3");

	        return make_tuple(newYaw, newPitch, newRoll);
	    }
	}
}


void JitterFromNetwork::setCurrentJitterStep(double endOfSimulation, double timeStep, double yaw, double pitch, double roll)
{

    // the old values are the new values from the last jitter step
    oldTimeStep = currentTimeStep;
    oldYaw = currentYaw;
    oldPitch = currentPitch;
    oldRoll = currentRoll;

    currentTimeStep = timeStep;
    currentYaw = deg2rad(yaw/3600.);
    currentPitch = deg2rad(pitch/3600.);
    currentRoll = deg2rad(roll/3600.);

	Log.info("JitterFromNetwork: jitter step set");
}


/**
 * \brief returns the jitter step (yaw, pitch, roll) at the given time value
 * 
 * \
 */
double JitterFromNetwork::getHeartbeatInterval()
{
	heartbeatInterval = currentTimeStep - oldTimeStep;

	return heartbeatInterval;
}
