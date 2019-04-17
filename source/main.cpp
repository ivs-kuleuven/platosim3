
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
    
    Log.info(string("PlatoSim ") + GIT_DESCRIBE);
    Log.info("Main: Log file includes 'error', 'warning', 'info', and 'debug'");

    std::mutex m;
    std::condition_variable condVar;

    // Initialise the simulation, and loop over all exposures using run()

    Simulation simulation(inputFilename, outputFilename, &m, &condVar);

    // check whether the simulation uses jitter from network

    if (simulation.getServerInstance() != NULL)
    {
	
    	Log.info("main: check");
    	
    	// the tcp connections have to run alongside the simulation so some threads have to be declared

    	std::thread serverThread(&TcpConnection::connectToServer, simulation.getServerInstance());
    	std::thread simulationThread(&Simulation::run, &simulation);

    	// gather the threads after completion and rejoin them

    	simulationThread.join();
    	serverThread.join();

    }
    else
    {
        simulation.run();
    }

    


    // That's it!

    logFile.close();
    return EXIT_SUCCESS;
}
