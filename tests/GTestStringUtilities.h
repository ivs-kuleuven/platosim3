#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "StringUtilities.h"

using namespace std;

TEST(StringUtilitiesTest, dtos) 
{

    using StringUtilities::dtos;

    LOG_STARTING_OF_TEST

    EXPECT_EQ("3.141593", dtos(3.1415926535897932385));
    EXPECT_EQ("3.142", dtos(3.1415926535897932385, false, 3));
    EXPECT_EQ("-3.141593", dtos(-3.1415926535897932385));

    EXPECT_EQ("3.141593e+00", dtos(3.1415926535897932385, true));
    EXPECT_EQ("-3.141593e+00", dtos(-3.1415926535897932385, true));

    EXPECT_EQ("2123456789.000000", dtos(2123456789.0));
    EXPECT_EQ("0.000002", dtos(0.000002123456789));  // I think we want to see 2.123456e-06 here

    EXPECT_EQ("2.123457e+09", dtos(2123456789.0, true));
    EXPECT_EQ("2.123457e-06", dtos(0.000002123456789, true));

    EXPECT_EQ("-2123456789.000000", dtos(-2123456789.0));
    EXPECT_EQ("-2123456789.00", dtos(-2123456789.0, false, 2));
    EXPECT_EQ("-0.000002", dtos(-0.000002123456789));  // I think we want to see -2.123456e-06 here

    EXPECT_EQ("-2.123457e+09", dtos(-2123456789.0, true));
    EXPECT_EQ("-2.123457e-06", dtos(-0.000002123456789, true));

    EXPECT_EQ("0.000000", dtos(0.0));
    EXPECT_EQ("0.000000e+00", dtos(0.0, true));
}







TEST(StringUtilitiesTest, environment)
{
    using StringUtilities::replaceEnvironmentVariable;

    LOG_STARTING_OF_TEST

    // Test what happends when no pattern is in the inputString

    string str = "Nothing here to replace";
    EXPECT_EQ("Nothing here to replace", replaceEnvironmentVariable(str));

    // Test what happens when the environment variable is not set/known

    str = "ENV['UNKNOWN_ENV']";
    EXPECT_EQ("ENV['UNKNOWN_ENV']", replaceEnvironmentVariable(str));

    // Make sure the environment variable is known to the tests

    setenv("PLATOSIM_PROJECT_HOME", "/Users/rik/Git/PlatoSim3", 1);

    // Check that the pattern is properly replaced

    str = "ENV['PLATOSIM_PROJECT_HOME']";
    EXPECT_EQ("/Users/rik/Git/PlatoSim3", replaceEnvironmentVariable(str));

    str = "ENV['PLATOSIM_PROJECT_HOME']/inputfiles";
    EXPECT_EQ("/Users/rik/Git/PlatoSim3/inputfiles", replaceEnvironmentVariable(str));

}




