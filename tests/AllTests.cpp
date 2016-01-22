
// Include all the tests here

#include "GTestConfigurationParameters.h"
#include "GTestFileUtilities.h"
#include "logger.h"

Logger Log;

int main(int argc, char **argv) 
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

