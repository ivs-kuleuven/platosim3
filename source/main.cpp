
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "Simulation.h"
#include "Logger.h"
#include "StringUtilities.h"
#include "version.h"

#include <thread>
#include <omp.h>


using namespace std;


Logger Log;


int main(int Narguments, char* arguments[])
{   
    using StringUtilities::ends_with;

    bool multipleSimulations;
    bool invalidInput = false;

    string firstArgument;

    if (Narguments >= 2) 
    {
        firstArgument = arguments[1];

        if (ends_with(arguments[1], "-version"))
        {
            cout << "PlatoSim " << GIT_DESCRIBE << endl;
            exit(EXIT_SUCCESS);
        }
        else if (firstArgument.substr(firstArgument.find_last_of(".") + 1) == "txt")
        {
            multipleSimulations = true;
        }
        else if (firstArgument.substr(firstArgument.find_last_of(".") + 1) == "yaml")
        {
            multipleSimulations = false;

            if (Narguments < 3 || (Narguments > 5))
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
        cerr << "Usage: platosim <inputfile> <outputfile> [<logfile>] [<logLevel>]" << endl;
        cerr << "       platosim -version" << endl;
        cerr << "       platosim <path_to_inputfile_list>" << endl;
        exit(EXIT_FAILURE); 
    }

    std::cout << invalidInput << std::endl;
    std::cout << multipleSimulations << std::endl;

    string inputFilename(arguments[1]);
    
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


    if (multipleSimulations == false)
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

        std::vector<Simulation*> simulationInstanceVec;

        // create Simulation objects as long as you have valid input- and output files
        for (int n = 0; n < inputFiles.size(); n++)
        {
            // create a new Simulation object
            Simulation* simulationInstance = new Simulation(inputFiles.at(n), outputFiles.at(n));

            // put the object handler in the respective vector
            simulationInstanceVec.emplace_back(simulationInstance);
        }

        std::vector<std::thread> simulationThreadVec;

        #pragma omp parallel for
        for (int n = 0; n < simulationInstanceVec.size(); n++)
        {        
            std::cout << "used threads: " << omp_get_num_threads() << std::endl;
    
            std::cout << "max threads: " << omp_get_max_threads() << std::endl;

            std::thread simulationThread(&Simulation::run, simulationInstanceVec.at(n));

            // put the thread handler in the respective vector
            simulationThreadVec.emplace_back(std::move(simulationThread));    
        }

        for (int n = 0; n < simulationThreadVec.size(); n++)
        {
            simulationThreadVec.at(n).join();
        }

        simulationInstanceVec.at(0)->deleteJitterAndDrift();

        //detach all Simulation Objects from the Jitter object and delete all class intances
        for (auto &i : simulationInstanceVec)
        {
           delete i;
        }
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
