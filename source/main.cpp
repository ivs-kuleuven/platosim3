
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
        if ((logLevel < 0) || (logLevel > 3)) 
        {
            cerr << "Error: logLevel should be either 0 (only errors), 2, or 3 (most verbose)" << endl;
            exit(EXIT_FAILURE);
        }
    }

    // Set up the log file


    ofstream logFile(logFilename);
    switch (logLevel)
    {
        case 0: Log.addOutputStream(logFile, ERROR);
                Log.addOutputStream(cerr,    ERROR);
                break;
        case 1: Log.addOutputStream(logFile, ERROR | WARNING);
                Log.addOutputStream(cerr,    ERROR | WARNING);
                break;
        case 2: Log.addOutputStream(logFile, ERROR | WARNING | INFO); 
                Log.addOutputStream(cerr,    ERROR | WARNING);                           // No excessive logging to stderr
                break;
        case 3: Log.addOutputStream(logFile, ERROR | WARNING | INFO | DEBUG);
                Log.addOutputStream(cerr,    ERROR | WARNING);                           // No excessive logging to stderr
                break;   
    }
    
    Log.info(string("PlatoSim ") + GIT_DESCRIBE);
    Log.info("Main: Log file includes 'error', 'warning', 'info', and 'debug'");


    // Initialise the simulation, and loop over all exposures using run()

    Simulation simulation(inputFilename, outputFilename);
    simulation.run();


    // That's it!

    logFile.close();
    return EXIT_SUCCESS;
}
