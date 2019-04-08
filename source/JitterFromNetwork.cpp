#include "JitterFromNetwork.h"



/**
 * \brief Constructor of the jitter from network class
 * 
 *       
 */
JitterFromNetwork::JitterFromNetwork(ConfigurationParameters &configParams)
{
    // Set the configuration parameters

    configure(configParams);

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