#include <fstream>

// Include all the tests here

#include "Logger.h"
#include "GTestConfigurationParameters.h"
#include "GTestFileUtilities.h"
#include "GTestStringUtilities.h"
#include "GTestPointSpreadFunction.h"
#include "GTestCamera.h"
#include "GTestDetector.h"
#include "GTestUnits.h"
#include "GTestArrayOperations.h"
#include "GTestPolynomial.h"
#include "GTestSkyCoordinates.h"
#include "PrettyPrinters.h"

Logger Log;


int main(int argc, char **argv) 
{
    ofstream logFile("log.txt");
    Log.addOutputStream(logFile, WARNING | ERROR | DEBUG | INFO);
    // Log.addOutputStream(cout, WARNING | ERROR | DEBUG | INFO);  // uncomment this line to show log messages on the Console
    Log.info("Main: Unit Test Log file includes 'warning', 'error', 'debug', and 'info'");

    ::testing::InitGoogleTest(&argc, argv);
    int returnValue = RUN_ALL_TESTS();

    logFile.close();
    return returnValue;
}

