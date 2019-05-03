
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
    
    // declare some variables for the inner thread communication
    std::mutex mServer;
    std::mutex mClient;
    std::condition_variable condVarServer;
    std::condition_variable condVarClient;

    // Initialise the simulation, and loop over all exposures using run()

    Simulation simulation(inputFilename, outputFilename, &mServer, &condVarServer, &mClient, &condVarClient);


    // declare and initialize the thread pointers
    std::thread* serverThread = NULL;
    std::thread* clientThread = NULL;
    std::thread* simulationThread =NULL;

    // declare a vector of pointers to the threads
    std::vector<std::thread*> threadVec;

    // depending on whether jitter from server or sending imagettes to client is active start the according threads

    // jitterFromNetwork AND sendImagettes is true
    if (simulation.getServerInstance() != NULL && simulation.getClientInstance() != NULL)
    {
        // create a server and a client thread

        serverThread = new std::thread(&TcpConnection::connectToServer, simulation.getServerInstance());
        threadVec.push_back(serverThread);

        clientThread = new std::thread(&TcpConnection::connectToClient, simulation.getClientInstance());
        threadVec.push_back(clientThread);
    }
    // jitterFromNetwork is true
    else if (simulation.getServerInstance() != NULL && simulation.getClientInstance() == NULL)
    {
        // create only the server thread

        serverThread = new std::thread(&TcpConnection::connectToServer, simulation.getServerInstance());
        threadVec.push_back(serverThread);
    }
    // sendImagettes is true
    else if (simulation.getServerInstance() == NULL && simulation.getClientInstance() != NULL)
    {
        // create only the client thread

        clientThread = new std::thread(&TcpConnection::connectToClient, simulation.getClientInstance());
        threadVec.push_back(clientThread);
    }

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
