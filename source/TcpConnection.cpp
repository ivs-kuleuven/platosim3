#include "TcpConnection.h"

/**
 * \brief Constructor: Constructor of the TcpConnection class
 * 
 *		  
 */

TcpConnection::TcpConnection(ConfigurationParameters &configParam, JitterGenerator* jitterFromNetwork)
{
	// declare a pointer to the jitterFromNetwork innstance

	jitterInstance = jitterFromNetwork;

	endOfSimulation = false;

}


/**
 * \brief Configure the TcpConnection object using the input parameter file
 * 
 * \param configParams  Contains all configuration parameters from the input file
 */

void TcpConnection::configure(ConfigurationParameters &configParams)
{

}

/**
 * \brief Destructor
 *  
 */

TcpConnection::~TcpConnection()
{
	
}



/**
 * \brief Function to be carried out in a thread parallel to the simulation thread. It connects to a server and gets jitter data from it.
 *  
 */
void TcpConnection::connectToServer()
{
	// connect to the server

	
	// repeat until the simulation is over

	while(!endOfSimulation)
	{

		// declare a lock


		// wait for a notification from the simulation


		// send a request to the server for the next jitter step

		
		// get the reply from server

		
		// process the jitter data


		// set the jitter step in the jitter generator object


		// unlock the simulation thread

		
		// check whether the simulation is at its end


	}

}



/**
 * \brief Function to be carried out in a thread parallel to the simulation thread. It connects to a client and send imagette data to it.
 *  
 */
void TcpConnection::connectToClient()
{

}
