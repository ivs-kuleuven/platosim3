#include "TcpConnection.h"

/**
 * \brief Constructor: Constructor of the TcpConnection class
 * 
 *		  
 */

TcpConnection::TcpConnection(ConfigurationParameters &configParam, std::condition_variable* cond_var, std::mutex* m, bool* notified, bool* newStep)
{
	// set some starting parameters

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
	Log.info("TcpConnection: jitter thread created");

	// declare socket

	zmq::context_t context(1);
	zmq::socket_t socket(context, ZMQ_REQ);

	// connect to the server

	socket.connect(tcpAddress);
	
	// repeat until the simulation is over

	while(!endOfSimulation)
	{

		//*newStepPointer = false;

		// declare a lock

		 std::unique_lock<std::mutex> lock(*mutexPointer);

		// wait for a notification from the simulation

		Log.info("TcpConnection: wait for notification from simulation thread");

		while(!*notifiedPointer)
        {
           	condVarPointer->wait(lock);
        }

        // if the thread is notified, a request to the server should be send

        // reset to the initial parameters
        *notifiedPointer = false;
        lock.unlock();
		bool validStep = false;

		Log.info("TcpConnection: request jitter step");
	
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
	
	           	Log.info("TcpConnection: jitter time step: " + to_string(currentJitterStepVec.at(1)));

	           	*newStepPointer = true;

	           	Log.info("TcpConnection: got new jitter step");
	
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


/**
 * \brief The server reply is in a string format - this function converts it to a vector of doubles
 * \	  TODO: build in some sanity checks
 *  
 */
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
 * \brief Function to be carried out in a thread parallel to the simulation thread. It connects to a client and sends imagette data to it.
 *  
 */
void TcpConnection::connectToClient()
{
	Log.info("TcpConnection: client thread created");

	// declare socket

	zmq::context_t context(1);
	zmq::socket_t socket(context, ZMQ_REQ);

	// connect to the client

	socket.connect(tcpAddressClient);

	// set the exposure counter
	int exposureCounter = 0;
	
	// repeat until the simulation is over

	while(!endOfSimulation)
	{
		// declare a lock

		std::unique_lock<std::mutex> lock(*mutexPointer);

		// wait for a notification from the simulation

		Log.info("TcpConnection: wait for new imagette notification from simulation thread");

		while(!*notifiedPointer)
       	{	
        	condVarPointer->wait(lock);
        }

        *notifiedPointer = false;
       	lock.unlock();

		// get the imagette from the detector object

		Log.info("TcpConnection: get imagette from detector");

		arma::Mat<float>* pixelMapPointer = detectorInstance->getCurrentPixelMap();

		exposureCounter++;

		if (exposureCounter == numExposures)
		{
			endOfSimulation = true;
		}

		// change the arma matrix to a more suitable format to send

		Log.info("TcpConnection: convert imagette from mat to char*");

		const char* imagetteString = convertMatrixToChar(pixelMapPointer, endOfSimulation);

		// send it to the client with a time stamp / imagette number
		
		Log.info("TcpConnection: send imagette string to client");

		zmq::message_t imagetteToSend (strlen(imagetteString));
	    memcpy(imagetteToSend.data(), imagetteString, strlen(imagetteString));
		socket.send(imagetteToSend);

		// notify the simulation thread to carry on the simulation

		Log.info("TcpConnection: notify simulation thread");

		*newStepPointer = true;

	   	condVarPointer->notify_one();

		// get the reply from client (this has to be done or zeroMQ won't work)
		zmq::message_t reply;
	    socket.recv(&reply);
	}
}

/**
 * \brief converts the endOfSimulation, rows, cols and the pixelmap values to a char (seperated by a blank space) to be send to the client
 *  
 */
const char* TcpConnection::convertMatrixToChar(arma::Mat<float>* pixelMapPointer, bool endOfSimulation)
{
	// declare whether the simulation is to end and get the rows and cols from the pixelMap

	int end = endOfSimulation;

	int rows = pixelMapPointer->n_rows;

	int cols = pixelMapPointer->n_cols;

	// write the values seperated by a white space to the string
	std::string imagetteString = to_string(end) + " " + to_string(rows) + " " + to_string(cols) + " ";

	// write every value of the pixelMap to the string
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			imagetteString = imagetteString + to_string(pixelMapPointer->at(i, j)) + " ";
		}
	}

	const char* imagetteChar = imagetteString.c_str();

	return imagetteChar;
}
