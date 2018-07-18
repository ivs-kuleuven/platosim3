#include "TcpConnection.h"


/**
 * \brief Constructor: This takes over the task of creating a jitter instance. since both are implemented
 * 		  with the singleton pattern, following tries to create an instance will get a pointer to those created here instead
 * 
 *		  
 */

TcpConnection::TcpConnection(string inputFilename)
{
	ConfigurationParameters configParams(inputFilename);

	configure(configParams);

	internalTime = -1.0;

	newJitterStep = false;

    endOfSimulation = true;

    client = false;

    // Depending on what the user requested, define the proper platform jitter generator

    if (!useJitter)
    {
        jitterGenerator = NoJitter::Instance();
    }
    else
    {
        if (jitterSource == "FromFile")
        {
            jitterGenerator = JitterFromFile::Instance(configParams);
        }
        else if (jitterSource == "FromRedNoise")
        {
            jitterGenerator = JitterFromRedNoise::Instance(configParams);
        }
        else if (jitterSource == "FromNetwork")
        {
            jitterGenerator = JitterFromNetwork::Instance(configParams, this);

            client = true;

            // this is to make sure, that the end of the simulation is determined by the server and not by platosim
            endOfSimulation = false;
        }
        else
        {
            string errorMessage = "Simulation: Jitter Source '" + jitterSource + "' is not supported.";
            Log.error(errorMessage);
            throw IllegalArgumentException(errorMessage);
        }
    }

    // whether or not, PlatoSim should send created imagettes to a client

    server = false; // TODO: implement sending imagettes 
}



/**
 * \brief Configure the Simulation object using the input parameter file
 * 
 * \param configParams  Contains all configuration parameters from the input file
 */

void TcpConnection::configure(ConfigurationParameters &configParams)
{
    useJitter         = configParams.getBoolean("Platform/UseJitter");
    jitterSource      = configParams.getString("Platform/JitterSource");
    tcpAddress        = configParams.getString("Platform/TcpAddress");

}


/**
 * \brief Destructor
 *  
 */

TcpConnection::~TcpConnection()
{
	
}


/**
* \brief function that runs in a seperate thread parallel to the simulator
* 		 which entertains an ongoing connection to a sender or receiver of data
*/
void TcpConnection::connectToServer(std::condition_variable* cond_var, bool* notified, bool* newStep, std::mutex* m)
{
    // only execute this thread, if platosim gets its jitter steps from a server
    if (client)
    {
        zmq::context_t context(1);
        zmq::socket_t socket(context, ZMQ_REQ);

        socket.connect(tcpAddress);

        while(!endOfSimulation)
        {
            // declare a lock 
            std::unique_lock<std::mutex> lock(*m);

            // wait for a notification
            while(!*notified)
            {
                cond_var->wait(lock);
            }

            // if the thread is notified, a request to the server should be send

            // reset to the initial parameters
            *notified = false;
            lock.unlock();

            // send a request to server
            zmq::message_t request (5);
            memcpy(request.data(), "New Jitter step, please.", 5);

            socket.send(request);

            // get the reply
            zmq::message_t reply;
            socket.recv(&reply);

            // fracture the received string
            std::string rpl = std::string(static_cast<char*>(reply.data()), reply.size());

            std::stringstream ss(rpl);

            double i;

            while (ss >> i)
            {
                currentJitterStepVec.push_back(i);

                if (ss.peek() == ',' || ss.peek() == ' ')
                {
                    ss.ignore();
                }
            }

            // check whether this is a new string
            if (currentJitterStepVec.at(1) > internalTime)
            {
                jitterGenerator->setCurrentJitterStep(currentJitterStepVec.at(0), currentJitterStepVec.at(1), currentJitterStepVec.at(2), currentJitterStepVec.at(3), currentJitterStepVec.at(4));

                *newStep = true;

                cond_var->notify_one();     
            }

            if (currentJitterStepVec.at(0) != 0)
            {
                endOfSimulation = true;
            }
            else
            {
                endOfSimulation = false;
            }

            currentJitterStepVec.clear();
        }
    }
}



/**
* \brief function that runs in a seperate thread parallel to the simulator
* 		 which entertains an ongoing connection to a sender or receiver of data
*/
void TcpConnection::connectToClient()
{
    // only execute this thread, if platosim gets its jitter steps from a server
    if (server)
    {
        // TODO: implement sending imagettes
    }
	
}


