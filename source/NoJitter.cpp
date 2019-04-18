#include "NoJitter.h"


/**
 * \brief function to acces the protected constructor, if there is no jitter instance yet
 */
JitterGenerator* NoJitter::Instance()
{
    if(_instance == 0)
    {
        _instance = new NoJitter();
    }
    return _instance;
}


/**
 * \brief Constructor
 * 
 */

NoJitter::NoJitter()
{
}










/**
 * \brief Destructor
 */

NoJitter::~NoJitter()
{

}







/**
 * \brief Configure this object using the parameters from the input parameters file
 * 
 * \param configParams  The configuration parameters
 */

void NoJitter::configure(ConfigurationParameters &configParams)
{
}









/**
 * \brief With No jitter, simply always return (yaw, pitch, roll)=(0,0,0)
 * \return (0, 0, 0)  [rad]
 */

tuple<double, double, double> NoJitter::getNextYawPitchRoll(double time)
{

    return make_tuple(0.0, 0.0, 0.0);
}








/**
 * \brief Return the heartbeat interval. For this null generator this is in principle
 *        infinity, in practice the largest double value.
 *
 * \details The heartbeat interval is the jitter time interval which is set to a fraction of 
 *          the jitter time scale. so that the changes in (yaw, pitch, roll) can still be 
 *          reliably tracked.
 *          
 * \return  heartbeatInterval [s]
 */

double NoJitter::getHeartbeatInterval()
{
    return numeric_limits<double>::max();
}

