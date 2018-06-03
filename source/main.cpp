
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <string>

#include "Simulation.h"
#include "Logger.h"
#include "StringUtilities.h"
#include "version.h"
#include "Clock.h"


#include "vector"

#include "JitterGenerator.h"

#include "DriftGenerator.h"

using namespace std;


Logger Log;


int main(int Narguments, char* arguments[])
{
    using StringUtilities::ends_with;

    bool paraSimulation;
    bool invalidInput = false;

    string firstArgument;


    // check for valid input arguments and derive which kind of simulation is to be started
    if(Narguments >= 2)
    {
        firstArgument = arguments[1];

        if (ends_with(firstArgument, "-version"))
        {
            cout << "PlatoSim " << GIT_DESCRIBE << endl;
            exit(EXIT_SUCCESS);
        }
        else if (firstArgument.substr(firstArgument.find_last_of(".") + 1) == "txt")
        {
            paraSimulation = true;
        }
        else if (firstArgument.substr(firstArgument.find_last_of(".") + 1) == "yaml")
        {
            paraSimulation = false;

            if (Narguments < 3)
            {
                invalidInput = true;

            }
        }
        else
        {
            invalidInput = true;
        }
    }
    else
    {
        invalidInput = true;    
    }


    if (invalidInput)
    {
        cerr << "Usage: platosim <inputfile> <outputfile> [<logfile>]" << endl;
        cerr << "       platosim -version" << endl;
        cerr << "       platosim <path_to_inputfile_list>" << endl;
        exit(EXIT_FAILURE); 
    }


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

    
    string inputFilename = firstArgument;

    // how platosim is called should determine whether it is used the "classic" way, or the parallel way

    if (paraSimulation == false)
    {
        string outputFilename(arguments[2]);

        // Initialise the simulation, and loop over all exposures using run()

        Simulation simulation(inputFilename, outputFilename);
        simulation.run();
    }
    else
    {
        // create a vector of the inputfiles used (its length indicates, how many simulations will be conducted simulataniously)
        std::vector<string> inputFiles;

        // a vector of outputfile names respectively
        std::vector<string> outputFiles;

        // check, whether the path to the inputfile yields a usable result
        ifstream inputFileList(inputFilename);
        if (inputFileList.is_open())
        {
            string line;
            while(getline(inputFileList, line))
            {
                // fill the inputFiles vector with the file names within the given folder
                inputFiles.emplace_back(line);

                // crop the ending from the filename
                size_t lastindex = line.find_last_of(".");
                string rawName = line.substr(0, lastindex);

                // create an outputfile name from the inputfilename and save it within the vector
                string outputName = rawName + ".hdf5";

                outputFiles.emplace_back(outputName);
            }
        }
        else
        {
            cerr << "<path_to_input_file> can't be opened" << endl;
            exit(EXIT_FAILURE);
        }

        // create a clock instance which dictates the cycle time for the simulations
        Clock* clockInstance = new Clock(inputFiles.at(0));

        std::vector<Simulation*> simulationInstanceVec;

        // create Simulation objects as long as you have valid input- and output files
        for (int n = 0; n < inputFiles.size(); n++)
        {
            // create a new Simulation object
            Simulation* simulationInstance = new Simulation(inputFiles.at(n), outputFiles.at(n));

            clockInstance->attach(simulationInstance);

            // put the object handler in the respective vector
            simulationInstanceVec.emplace_back(simulationInstance);
        }

        std::cout << "Start Simulation" << std::endl;

        // start the Simulation
        clockInstance->startSimulation();

        std::cout << "End Simulation" << std::endl;

        //detach all Simulation Objects from the Jitter object and delete all class intances
        for (auto &i : simulationInstanceVec)
        {
           clockInstance->detach(i);
           delete i;
        }

        delete clockInstance;
    }

    // That's it!

    logFile.close();
    return EXIT_SUCCESS;
}


/**
 * \brief lazy initialization of the jitter instance
 */
JitterGenerator* JitterGenerator::_instance = 0;

 /**
 * \brief lazy initialization of the drift instance
 */
DriftGenerator* DriftGenerator::_instance = 0;
