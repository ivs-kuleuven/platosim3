
// Include all the tests here

#include "Logger.h"

#include "GTestConfigurationParameters.h"
#include "GTestFileUtilities.h"
#include "GTestPointSpreadFunction.h"

Logger Log;

int main(int argc, char **argv) 
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

