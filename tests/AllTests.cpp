#include <fstream>

// Include all the tests here

#include "Logger.h"

#include "GTestConfigurationParameters.h"
#include "GTestFileUtilities.h"
#include "GTestPointSpreadFunction.h"
#include "GTestPsfConvolution.h"

Logger Log;

int main(int argc, char **argv) 
{
    ofstream logFile("log.txt");
    Log.addOutputStream(logFile, WARNING | ERROR | DEBUG | INFO);
    Log.info("Main: Log file includes 'warning', 'error', 'debug', and 'info'");

    ::testing::InitGoogleTest(&argc, argv);
    int returnValue = RUN_ALL_TESTS();

    logFile.close();
    return returnValue;
}

