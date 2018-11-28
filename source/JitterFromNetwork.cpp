#include "JitterFromNetwork.h"


/**
 * \brief function to acces the protected constructor, if there is no jitter instance yet
 */
JitterGenerator* JitterFromNetwork::Instance(ConfigurationParameters &configParams, double readoutTimeBeforeNextExposure, TcpConnection* tcpConnection)
{
    if(_instance == 0)
    {
        _instance = new JitterFromNetwork(configParams, readoutTimeBeforeNextExposure, tcpConnection);
    }

    return _instance;
}


/**
 * \brief this is a variant of the intance function called by the simulation instances. this is, so the tcp connection object
 * \      does not need to be known by the simulation class
 */
JitterGenerator* JitterFromNetwork::Instance(ConfigurationParameters &configParams, double readoutTimeBeforeNextExposure)
{
    return _instance;
}


/**
 * \brief Constructor of the jitter from network class
 * 
 * \      
 */
JitterFromNetwork::JitterFromNetwork(ConfigurationParameters &configParams, double readoutTimeBeforeNextExposure, TcpConnection* tcpConnection)
: lastYaw(0.0), lastPitch(0.0), lastRoll(0.0), internalTime(0.0)
{
    // Set the configuration parameters

    configure(configParams);

    tcpConnectionInstance = tcpConnection;

    jitterTimeInterval = 0.0;

    oldTimeStep = 0.0;

    currentTimeStep = 0.0;

    oldYaw = 0.0;
    oldPitch = 0.0;
    oldRoll = 0.0;

    currentYaw = 0.0;
    currentPitch = 0.0;
    currentRoll = 0.0;

    endSimulation = false;
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
    // add parameters here later, if they're needed
}


/**
 * \brief Get the next (yaw, pitch, roll) values from a server
 * 
 * \note Also during CCD readout, the spacecraft jitters, to the user needs to take this into
 *       account when passing 'timeInterval'.
 * 
 * \param time   The time point for which the (yaw, pitch, roll) is requested [s]
 * 
 * \return (newYaw, newPitch, newRoll)  [rad]
 */

tuple<double, double, double> JitterFromNetwork::getNextYawPitchRoll(double time)
{

    if (time == currentTimeStep)
    {
        return make_tuple(currentYaw, currentPitch, currentRoll);
    }
    else
    {
        // interpolate the jitterstep dependent on the old and new jitterstep and the currentTime

        const double weight1 = (time - oldTimeStep) / (currentTimeStep - oldTimeStep);
        const double weight2 = (currentTimeStep - time) / (currentTimeStep - oldTimeStep);
        const double newYaw   = oldYaw   * weight2 + currentYaw   * weight1;
        const double newPitch = oldPitch * weight2 + currentPitch * weight1;
        const double newRoll  = oldRoll  * weight2 + currentRoll  * weight1;

        return make_tuple(newYaw, newPitch, newRoll);
    }
}


/**
 * \brief Return the heartbeat interval of this network jitter generator
 * 
 * \details The heartbeat interval is the difference between the time values of the newest jitterstep
 *          and the one before          
 *          
 * \return  heartbeatInterval [s]
 */

double JitterFromNetwork::getHeartbeatInterval()
{
    jitterTimeInterval = currentTimeStep - oldTimeStep;
    
    return jitterTimeInterval;
}

/**
 * \brief this function is called by the tcp connection client object, that sets the current jitter parameters, when it gets them 
 */

void JitterFromNetwork::setCurrentJitterStep(double endOfSimulation, double timeStep, double yaw, double pitch, double roll)
{
    if (endOfSimulation != 0)
    {
        endSimulation = true;
    }

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
 * \brief this function is called by the clock instance in a parallel simulation to determine, whether the simulation relies on jittersteps from a server  
 */
bool JitterFromNetwork::isClient()
{
    return true;
}

/**
 * \brief this function is called by the clock instance in a parallel simulation to determine, whether the server is still sending jittersteps  
 */

bool JitterFromNetwork::simulationEnd()
{
    return endSimulation;
}