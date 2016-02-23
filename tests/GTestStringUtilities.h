#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "StringUtilities.h"

using namespace std;

TEST(StringUtilitiesTest, dtos) 
{

    LOG_STARTING_OF_TEST

    EXPECT_EQ("3.141593", StringUtilities::dtos(3.1415926535897932385));
    EXPECT_EQ("3.142", StringUtilities::dtos(3.1415926535897932385, false, 3));
    EXPECT_EQ("-3.141593", StringUtilities::dtos(-3.1415926535897932385));

    EXPECT_EQ("3.141593e+00", StringUtilities::dtos(3.1415926535897932385, true));
    EXPECT_EQ("-3.141593e+00", StringUtilities::dtos(-3.1415926535897932385, true));

    EXPECT_EQ("2123456789.000000", StringUtilities::dtos(2123456789.0));
    EXPECT_EQ("0.000002", StringUtilities::dtos(0.000002123456789));  // I think we want to see 2.123456e-06 here

    EXPECT_EQ("2.123457e+09", StringUtilities::dtos(2123456789.0, true));
    EXPECT_EQ("2.123457e-06", StringUtilities::dtos(0.000002123456789, true));

    EXPECT_EQ("-2123456789.000000", StringUtilities::dtos(-2123456789.0));
    EXPECT_EQ("-2123456789.00", StringUtilities::dtos(-2123456789.0, false, 2));
    EXPECT_EQ("-0.000002", StringUtilities::dtos(-0.000002123456789));  // I think we want to see -2.123456e-06 here

    EXPECT_EQ("-2.123457e+09", StringUtilities::dtos(-2123456789.0, true));
    EXPECT_EQ("-2.123457e-06", StringUtilities::dtos(-0.000002123456789, true));

    EXPECT_EQ("0.000000", StringUtilities::dtos(0.0));
    EXPECT_EQ("0.000000e+00", StringUtilities::dtos(0.0, true));
}



