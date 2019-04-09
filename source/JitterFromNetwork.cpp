#include "JitterFromNetwork.h"



/**
 * \brief Constructor of the jitter from network class
 * 
 *       
 */
JitterFromNetwork::JitterFromNetwork(ConfigurationParameters &configParams, double readoutTimeBeforeNextExposure, std::condition_variable* cond_var, std::mutex* m, bool* notified, bool* newStep)
{
    // Set the configuration parameters

    configure(configParams);

    condVarPointer = cond_var;
    mutexPointer = m;
    notifiedPointer = notified;
    newStepPointer = newStep;

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
   
}



/**
 * \brief returns the jitter step (yaw, pitch, roll) at the given time value
 * 
 * \
 */
tuple<double, double, double> JitterFromNetwork::getNextYawPitchRoll(double time)
{
	const double yaw = 0.0;
	const double pitch = 0.0;
	const double roll = 0.0;

	return make_tuple(yaw, pitch, roll);
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
}



/**
 * \brief returns the jitter step (yaw, pitch, roll) at the given time value
 * 
 * \
 */
double JitterFromNetwork::getHeartbeatInterval()
{
	double hearbeatInterval = 0.0;
	return heartbeatInterval;
}
