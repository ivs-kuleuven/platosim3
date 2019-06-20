#include "TcpConnection.h"

/**
 * \brief Constructor: Constructor of the TcpConnection class
 * 
 *        
 */

TcpConnection::TcpConnection(ConfigurationParameters &configParam, std::condition_variable* cond_var, std::mutex* m, bool* notified, bool* newStep, bool* endSimulation)
{
    // set some starting parameters

    configure(configParam);

    condVarPointer = cond_var;
    mutexPointer = m;
    notifiedPointer = notified;
    newStepPointer = newStep;

    endOfSimulation = endSimulation;
}


/**
 * \brief Configure the TcpConnection object using the input parameter file
 * 
 * \param configParams  Contains all configuration parameters from the input file
 */

void TcpConnection::configure(ConfigurationParameters &configParams)
{
    internalTime = -1.0;

    tcpAddressInputServer = configParams.getString("ControlTcpConnection/TcpAddressInputServer");
    tcpAddressJitterServer = configParams.getString("ControlTcpConnection/TcpAddressJitterServer");
    tcpAddressImagetteClient = configParams.getString("ControlTcpConnection/TcpAddressImagetteClient");
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
void TcpConnection::connectToJitterServer(bool active)
{
    if (active)
    {
        Log.info("TcpConnection: jitter server thread created");

        // declare socket
        zmq::context_t context(1);
        zmq::socket_t socket(context, ZMQ_REQ);

        // connect to the server
        socket.connect(tcpAddressJitterServer);
        
        // repeat until the simulation is over
        while(!*endOfSimulation)
        {
            // declare a lock
            std::unique_lock<std::mutex> lock(*mutexPointer);

            // wait for a notification from the simulation
            Log.info("TcpConnection: wait for new jitter step request notification from simulation thread");

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
        
                    // notify the simulation thread that a new step is available
                    *newStepPointer = true;
        
                    Log.info("TcpConnection: got new jitter step");
            
                    condVarPointer->notify_one();
        
                    validStep = true;     
                }
                
                // check whether the simulation is at its end
                if (currentJitterStepVec.at(0) != 0)
                {
                   *endOfSimulation = true;
                }   
            }
        }

        Log.info("TcpConnection: get jitter thread ends");
    }
}


/**
 * \brief The server reply is in a string format - this function converts it to a vector of doubles
 * \      TODO: build in some sanity checks
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
void TcpConnection::connectToImagetteClient(bool active)
{
    if (active)
    {
        Log.info("TcpConnection: imagette client thread created");

        // declare socket
        zmq::context_t context(1);
        zmq::socket_t socket(context, ZMQ_REQ);

        // connect to the client
        socket.connect(tcpAddressImagetteClient);

        // repeat until the simulation is over
        bool endSimulation = false;

        while(!endSimulation)
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
            
            // if the thread is notified by the detector and the end of the simulation has been declared (by jitterFromNetwork or the Simulation)
            // this ends the loop and closes the thread as soon as the last imagette is sent
            if (*endOfSimulation)
            {
                Log.info("TcpConnection: Simulation end");
                endSimulation = true;
            }

            // get the imagette from the detector object
            Log.info("TcpConnection: get imagette from detector");

            arma::Mat<float>* pixelMapPointer = detectorInstance->getCurrentPixelMap();

            // change the arma matrix to a more suitable format to send
            Log.info("TcpConnection: convert imagette from mat to char*");

            const char* imagetteString = convertMatrixToChar(pixelMapPointer);

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

        Log.info("TcpConnection: send imagette thread ends");
    }
}

/**
 * \brief converts the endOfSimulation, rows, cols and the pixelmap values to a char (seperated by a blank space) to be send to the client
 *  
 */
const char* TcpConnection::convertMatrixToChar(arma::Mat<float>* pixelMapPointer)
{
    // declare whether the simulation is to end and get the rows and cols from the pixelMap
    int end = *endOfSimulation;

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





/**
 * \brief Function to be carried out in a thread parallel to the simulation thread. It connects to a server and gets input data from it.
 *  
 */
void TcpConnection::connectToInputServer(bool active)
{
    if (active)
    {
        Log.info("TcpConnection: input server thread created");

        // declare socket
        zmq::context_t context(1);
        zmq::socket_t socket(context, ZMQ_REQ);

        // connect to the server
        socket.connect(tcpAddressInputServer);
        
        // repeat until the simulation is over
        int numColumns = 0;
        int numRows = 0;

        bool endSimulation = false;

        while(!endSimulation)
        {
            // declare a lock
            std::unique_lock<std::mutex> lock(*mutexPointer);

            // wait for a notification from the simulation
            Log.info("TcpConnection: wait for input request notification from detector thread");

            while(!*notifiedPointer)
            {
                condVarPointer->wait(lock);
            }

            // if the thread is notified by the detector and the end of the simulation has been declared (by jitterFromNetwork or the Simulation)
            // this ends the loop and closes the thread as soon as the last imagette is sent
            if (*endOfSimulation)
            {
                Log.info("TcpConnection: Simulation end");
                endSimulation = true;
            }

            // reset to the initial parameters
            *notifiedPointer = false;
            lock.unlock();
            bool validStep = false;

            int imagetteNumber = detectorInstance->getCurrentImagetteNumber();

            // create a string composed of whether the simulation is at its end the current exposure number

            const char* imagetteNumberString = (std::to_string(int(endSimulation)) + " " + std::to_string(imagetteNumber)).c_str();

            std::cout << imagetteNumberString << std::endl;

            // send it to the client with a imagette number
            
            Log.info("TcpConnection: send imagette number string to client");

            zmq::message_t imagetteNumberToSend (strlen(imagetteNumberString));

            memcpy(imagetteNumberToSend.data(), imagetteNumberString, strlen(imagetteNumberString));

            socket.send(imagetteNumberToSend);

            // get the reply from server
            zmq::message_t reply;

            socket.recv(&reply);

            std::string replyString = std::string(static_cast<char*>(reply.data()), reply.size());

            std::vector<double> currentInputVec = processServerReply(replyString);

            // if the size of the subfield has changed, change the size and position values of the subfield in the detector 
            if (currentInputVec.at(0) != numColumns && currentInputVec.at(1) != numRows && !endSimulation)
            {
                detectorInstance->setInputParametersFromServer(currentInputVec.at(0), currentInputVec.at(1), currentInputVec.at(2), currentInputVec.at(3), currentInputVec.at(4), currentInputVec.at(5), currentInputVec.at(6));
            
                Log.info("TcpConnection: set input Parameters");
            }

            numColumns = currentInputVec.at(0);

            numRows = currentInputVec.at(1);


            Log.info("TcpConnection: got new input step");
            
            // notify the simulation thread
            *newStepPointer = true;
            
            condVarPointer->notify_one();

        }

        Log.info("TcpConnection: end input server thread");
    }
}


