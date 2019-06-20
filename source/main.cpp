
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "Simulation.h"
#include "Logger.h"
#include "StringUtilities.h"
#include "version.h"

using namespace std;


Logger Log;


int main(int Narguments, char* arguments[])
{
    using StringUtilities::ends_with;

    if (Narguments == 2 && ends_with(arguments[1], "-version"))
    {
        cout << "PlatoSim " << GIT_DESCRIBE << endl;
        exit(EXIT_SUCCESS);
    }

    // Platosim expects the filename of the configuration parameters.
    // Exit if this filename is not given.

    if ((Narguments < 3) || (Narguments > 5))
    {
        cerr << "Usage: platosim <inputfile> <outputfile> [<logfile>] [<logLevel>]" << endl;
        cerr << "       platosim -version" << endl;
        exit(EXIT_FAILURE);
    }

    string inputFilename(arguments[1]);
    string outputFilename(arguments[2]);
    string logFilename = "log.txt";
    if (Narguments >= 4)
    {
        logFilename = arguments[3];
    }

    int logLevel = 3;
    if (Narguments == 5)
    {
        logLevel = atoi(arguments[4]);
        if ((logLevel < 1) || (logLevel > 3)) 
        {
            cerr << "Error: logLevel should be either 1 (least verbose), 2, or 3 (most verbose)" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // Set up the log file

    Log.addOutputStream(cerr, ERROR | WARNING);

    ofstream logFile(logFilename);
    switch (logLevel)
    {
        case 1: Log.addOutputStream(logFile, ERROR | WARNING);
                break;
        case 2: Log.addOutputStream(logFile, ERROR | WARNING | INFO); 
                break;
        case 3: Log.addOutputStream(logFile, ERROR | WARNING | INFO | DEBUG);
                break;   
    }
    
    // declare some variables for the communication between threads 

    std::mutex mJitterServer;
    std::mutex mImagetteClient;
    std::mutex mInputServer;
    std::condition_variable condVarJitterServer;
    std::condition_variable condVarImagetteClient;
    std::condition_variable condVarInputServer;

    std::vector<std::tuple<std::mutex*, std::condition_variable*> > threadCommunicationVec;

    threadCommunicationVec = { std::make_tuple( &mJitterServer, &condVarJitterServer), std::make_tuple( &mImagetteClient, &condVarImagetteClient), std::make_tuple( &mInputServer, &condVarInputServer) };

    // Initialise the simulation object

    Simulation simulation(inputFilename, outputFilename, threadCommunicationVec);

    // declare and initialize the thread pointers
    std::thread* jitterServerThread = NULL;
    std::thread* imagetteClientThread = NULL;
    std::thread* inputServerThread = NULL;
    std::thread* simulationThread = NULL;

    // declare a vector of pointers to the threads
    std::vector<std::thread*> threadVec;

    // declare the attributes for the different thread methods
    bool jitterServerActive = false;
    bool inputServerActive = false;
    bool imagetteClientActive = false;

    // check whether the simulation is done with jitter from a server via tcp connection 
    if (simulation.getJitterServerInstance() != NULL)
    {
        jitterServerActive = true;
    }

    // check whether the simulation is done with input from a server via tcp connection
    if (simulation.getInputServerInstance() != NULL)
    {
        inputServerActive = true;
    }

    // check whether the simulation sends the created imagettes to a client via a tcp connection
    if (simulation.getImagetteClientInstance() != NULL)
    {
        imagetteClientActive = true;
    }

    // create the jitter server thread
    jitterServerThread = new std::thread(&TcpConnection::connectToJitterServer, simulation.getJitterServerInstance(), jitterServerActive);
    threadVec.push_back(jitterServerThread);

    // create the input server thread
    inputServerThread = new std::thread(&TcpConnection::connectToInputServer, simulation.getInputServerInstance(), inputServerActive);
    threadVec.push_back(inputServerThread);

    // create the imagette client thread
    imagetteClientThread = new std::thread(&TcpConnection::connectToImagetteClient, simulation.getImagetteClientInstance(), imagetteClientActive);
    threadVec.push_back(imagetteClientThread);

    // create the simulation thread
    simulationThread = new std::thread(&Simulation::run, &simulation);
    threadVec.push_back(simulationThread);

    // after the simulation join all threads and delete them
    for (auto &i : threadVec)
    {
        i->join();
        delete(i);
    }

    // That's it!

    logFile.close();

    return EXIT_SUCCESS;
}
