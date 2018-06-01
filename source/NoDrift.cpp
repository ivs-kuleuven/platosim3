#include "NoDrift.h"


/**
 * \brief function to acces the protected constructor, if there is no drift instance yet
 */
DriftGenerator* NoDrift::Instance()
{
    if(_instance == 0)
    {
        _instance = new NoDrift();
    }
    return _instance;
}




/**
 * \brief Constructor
 * 
 * \param configParams The configuration parameters from the input parameters file
 */

NoDrift::NoDrift()
{
}










/**
 * \brief Destructor
 */

NoDrift::~NoDrift()
{

}







/**
 * \brief Configure this object using the parameters from the input parameters file
 * 
 * \param configParams  The configuration parameters
 */

void NoDrift::configure()
{
}









/**
 * \brief With No Drift, simply always return (yaw, pitch, roll)=(0,0,0)
 * \return (0, 0, 0)  [rad]
 */

tuple<double, double, double> NoDrift::getNextYawPitchRoll(double time)
{
    return make_tuple(0.0, 0.0, 0.0);
}








/**
 * \brief Return the heartbeat interval. For this null generator this is in principle
 *        infinity, in practice the largest double value.
 * 
 * \details The heartbeat interval is the drift time interval which is set to a fraction of 
 *          the drift time scale. so that the changes in (yaw, pitch, roll) can still be 
 *          reliably tracked.
 *          
 * \return  heartbeatInterval [s]
 */

double NoDrift::getHeartbeatInterval()
{
    return numeric_limits<double>::max();
}

