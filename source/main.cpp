
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "Simulation.h"
#include "Logger.h"


using namespace std;


Logger Log;


int main(int Narguments, char* arguments[])
{
    // Platosim expects the filename of the configuration parameters.
    // Exit if this filename is not given.

    if (Narguments != 3)
    {
        cerr << "Usage: platosim <inputfile> <outputfile>" << endl;
        exit(EXIT_FAILURE);
    }

    string inputFilename(arguments[1]);
    string outputFilename(arguments[2]);


    // Set up the log file

    ofstream logFile("log.txt");
    Log.addOutputStream(cerr,    WARNING | ERROR);
    Log.addOutputStream(logFile, WARNING | ERROR | DEBUG | INFO);
    Log.info("Main: Log file includes 'warning', 'error', 'debug', and 'info'");


    // Initialise the simulation, and loop over all exposures using run()

    Simulation simulation(inputFilename, outputFilename);
    simulation.run();


    // That's it!

    logFile.close();
    return EXIT_SUCCESS;
}