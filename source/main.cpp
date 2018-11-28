
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <string>

#include <condition_variable>

#include "Simulation.h"
#include "Logger.h"
#include "StringUtilities.h"
#include "version.h"
#include "Clock.h"


#include "vector"

#include "JitterGenerator.h"

#include "DriftGenerator.h"

#include "TcpConnection.h"

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

    
    string inputFilename = firstArgument;

    // how platosim is called should determine whether it is used the "classic" way, or the parallel way

    if (paraSimulation == false)
    {
        string outputFilename(arguments[2]);

        // Initialise the simulation, and loop over all exposures using run()

        Simulation simulation(inputFilename, outputFilename, paraSimulation);

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

        // create tcp connection instances for a server and a client object
        TcpConnection* serverInstance = new TcpConnection(inputFiles.at(0));

        TcpConnection* clientInstance = new TcpConnection(inputFiles.at(0));


        std::vector<Simulation*> simulationInstanceVec;

        // create Simulation objects as long as you have valid input- and output files
        for (int n = 0; n < inputFiles.size(); n++)
        {
            // create a new Simulation object
            Simulation* simulationInstance = new Simulation(inputFiles.at(n), outputFiles.at(n), paraSimulation);

            // attach the simulation instance as observer to the clock instance
            clockInstance->attach(simulationInstance);

            // put the object handler in the respective vector
            simulationInstanceVec.emplace_back(simulationInstance);
        }

        // create a conditional variable object instance
        std::condition_variable cond_var;
        std::condition_variable* conPtr = &cond_var;

        bool notified = false;
        bool* pNotified = &notified;

        bool newStep = false;
        bool* pNewStep = &newStep;

        std::mutex m;
        std::mutex* pM = &m;

        std::thread simulationThread(&Clock::startSimulation, clockInstance, conPtr, pNotified, pNewStep, pM);

        std::thread serverThread(&TcpConnection::connectToClient, serverInstance);

        std::thread clientThread(&TcpConnection::connectToServer, clientInstance, simulationInstanceVec.at(0), conPtr, pNotified, pNewStep, pM);

        serverThread.join();

        clientThread.join(); 

        simulationThread.join();

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
 * \brief needed for the lazy initialization of the jitter instance
 */
JitterGenerator* JitterGenerator::_instance = 0;

 /**
 * \brief needed for the lazy initialization of the drift instance
 */
DriftGenerator* DriftGenerator::_instance = 0;
