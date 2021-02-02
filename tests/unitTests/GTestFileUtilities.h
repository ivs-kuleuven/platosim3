#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "FileUtilities.h"

using namespace std;


TEST(FileUtilitiesTest, fileExists) 
{
    
    LOG_STARTING_OF_TEST

    char * env = getenv("PLATO_PROJECT_HOME");
    ASSERT_FALSE(env == NULL);

    ASSERT_TRUE(FileUtilities::fileExists(string(getenv("PLATO_PROJECT_HOME")) + "/testData/input_ConfigurationParametersTest.yaml"));

    ASSERT_FALSE(FileUtilities::fileExists(string(getenv("PLATO_PROJECT_HOME")) + "/wrongDirectory/input.yaml"));

}


TEST(FileUtilitiesTest, isRelative) 
{
    
    LOG_STARTING_OF_TEST

    char * env = getenv("PLATO_PROJECT_HOME");
    ASSERT_FALSE(env == NULL);

    ASSERT_FALSE(FileUtilities::isRelative(string(getenv("PLATO_PROJECT_HOME")) + "/testData/input_ConfigurationParametersTest.yaml"));

    ASSERT_FALSE(FileUtilities::isRelative("/Users/rik/Git/PlatoSim3/testData/input.yaml"));

}

