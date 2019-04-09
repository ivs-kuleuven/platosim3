#include "TcpConnection.h"

/**
 * \brief Constructor: Constructor of the TcpConnection class
 * 
 *		  
 */

TcpConnection::TcpConnection(ConfigurationParameters &configParam, JitterGenerator* jitterFromNetwork, std::condition_variable* cond_var, std::mutex* m, bool* notified, bool* newStep)
{
	// set some starting parameters

	jitterInstance = jitterFromNetwork;

	configure(configParam);

	condVarPointer = cond_var;
	mutexPointer = m;
	notifiedPointer = notified;
	newStepPointer = newStep;
}


/**
 * \brief Configure the TcpConnection object using the input parameter file
 * 
 * \param configParams  Contains all configuration parameters from the input file
 */

void TcpConnection::configure(ConfigurationParameters &configParams)
{
	endOfSimulation = false;

	internalTime = -1.0;

	tcpAddress = configParams.getString("Platform/TcpAddress");
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
	// declare socket

	zmq::context_t context(1);
	zmq::socket_t socket(context, ZMQ_REQ);

	// connect to the server

	socket.connect(tcpAddress);
	
	// repeat until the simulation is over

	while(!endOfSimulation)
	{

		// declare a lock

		 std::unique_lock<std::mutex> lock(*mutexPointer);

		// wait for a notification from the simulation

		while(!*notifiedPointer)
            	{
                	condVarPointer->wait(lock);
            	}

                // if the thread is notified, a request to the server should be send

               	// reset to the initial parameters
               	*notifiedPointer = false;
               	lock.unlock();
		bool validStep = false;

		// request new jitter steps until a valid step is send
		while(!validStep)
		{
			// send a request to the server for the next jitter step
	
	                zmq::message_t request (5);
	                memcpy(request.data(), "New Jitter step, please.", 5);
			socket.send(request);
			
			// get the reply from server
			zmq::message_t reply;
	                socket.recv(&reply);
			
			// process the jitter data

			std::string replyString = std::string(static_cast<char*>(reply.data()), reply.size());
	
			std::vector<double> currentJitterStepVec = processServerReply(replyString);
	

			// set the jitter step in the jitter generator object if its time stamp is higher than the internal time
	
	            	if (currentJitterStepVec.at(1) > internalTime)
	            	{
	                	jitterInstance->setCurrentJitterStep(currentJitterStepVec.at(0), currentJitterStepVec.at(1), currentJitterStepVec.at(2), currentJitterStepVec.at(3), currentJitterStepVec.at(4));
	
	                	*newStepPointer = true;
	
	                	condVarPointer->notify_one();

				validStep = true;     
	            	}
			
			// check whether the simulation is at its end
	
			if (currentJitterStepVec.at(0) != 0)
	           	{
	                	endOfSimulation = true;
	            	}	
	            	else
	            	{
	                	endOfSimulation = false;
            		}
		}
	}
}

std::vector <double> TcpConnection::processServerReply(string replyString)
{
	// fracture the received string

        std::stringstream ss(replyString);

	std::vector<double> currentJitterStepVec;

        double i;

	// convert the string to double values

        while (ss >> i)
        {
             currentJitterStepVec.push_back(i);

             if (ss.peek() == ',' || ss.peek() == ' ')
             {
                 ss.ignore();
             }
         }
	
	return currentJitterStepVec;

}




/**
 * \brief Function to be carried out in a thread parallel to the simulation thread. It connects to a client and send imagette data to it.
 *  
 */
void TcpConnection::connectToClient()
{

}
