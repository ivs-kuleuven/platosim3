#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "FileUtilities.h"

using namespace std;


TEST(FileUtilitiesTest, fileExists) 
{
    
    LOG_STARTING_OF_TEST

    ASSERT_TRUE(FileUtilities::fileExists("../testData/input_ConfigurationParametersTest.yaml"));

    ASSERT_FALSE(FileUtilities::fileExists("../wrongDirectory/input.yaml"));

}


TEST(FileUtilitiesTest, isRelative) 
{
    
    LOG_STARTING_OF_TEST

    ASSERT_TRUE(FileUtilities::isRelative("../testData/input_ConfigurationParametersTest.yaml"));

    ASSERT_FALSE(FileUtilities::isRelative("/Users/rik/Git/PlatoSim3/testData/input.yaml"));

}

