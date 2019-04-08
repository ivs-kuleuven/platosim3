#include "TcpConnection.h"

/**
 * \brief Constructor: Constructor of the TcpConnection class
 * 
 *		  
 */

TcpConnection::TcpConnection(string inputFilename)
{


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

}



/**
 * \brief Function to be carried out in a thread parallel to the simulation thread. It connects to a client and send imagette data to it.
 *  
 */
void TcpConnection::connectToClient()
{

}
