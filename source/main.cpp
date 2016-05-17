
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

    if (Narguments < 3)
    {
        cerr << "Usage: platosim <inputfile> <outputfile> [<logfile>]" << endl;
        cerr << "       platosim -version" << endl;
        exit(EXIT_FAILURE);
    }

    string inputFilename(arguments[1]);
    string outputFilename(arguments[2]);
    string logFilename = "log.txt";
    if (Narguments == 4)
    {
        logFilename = arguments[3];
    }

    // Set up the log file

    ofstream logFile(logFilename);
    Log.addOutputStream(cerr,    WARNING | ERROR);
    Log.addOutputStream(logFile, WARNING | ERROR | DEBUG | INFO);
    Log.info(string("PlatoSim ") + GIT_DESCRIBE);
    Log.info("Main: Log file includes 'warning', 'error', 'debug', and 'info'");


    // Initialise the simulation, and loop over all exposures using run()

    Simulation simulation(inputFilename, outputFilename);
    simulation.run();


    // That's it!

    logFile.close();
    return EXIT_SUCCESS;
}