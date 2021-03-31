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


TEST(StringUtilitiesTest, replaceEnvironmentVariable)
{
    using StringUtilities::replaceEnvironmentVariable;

    LOG_STARTING_OF_TEST

    // if no ENV['var'] pattern, the string should just be returned

    EXPECT_EQ("ABC__XYZ__DEF", replaceEnvironmentVariable("ABC__XYZ__DEF"));

    // no environment variable __XYZ__ should exist at this point

    EXPECT_EQ("ABCENV['__XYZ__']DEF", replaceEnvironmentVariable("ABCENV['__XYZ__']DEF"));

    // Make sure the environment variable __XYZ__ is known to the tests

    setenv("__XYZ__", "+++QWERTY+++", 1);

    string str = "ENV['__XYZ__']";
    EXPECT_EQ("+++QWERTY+++", replaceEnvironmentVariable(str));

    str = "ABCENV['__XYZ__']DEF";
    EXPECT_EQ("ABC+++QWERTY+++DEF", replaceEnvironmentVariable(str));

    str = "ENVENV['__XYZ__']ENV";
    EXPECT_EQ("ENV+++QWERTY+++ENV", replaceEnvironmentVariable(str));

    // Test the condition that ENV['var'] pattern is not complete (e.g. typo!)
 
    str = "ENVENV['__XYZ__'}";
    EXPECT_EQ(str, replaceEnvironmentVariable(str));

}




