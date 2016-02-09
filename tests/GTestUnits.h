#include "gtest/gtest.h"

#include "Units.h"
#include "Constants.h"

using namespace std;

TEST(UnitsTest, deg2rad2deg)
{
    EXPECT_DOUBLE_EQ(180. / Constants::PI, rad2deg(1.0));
    EXPECT_DOUBLE_EQ(Constants::PI / 180., deg2rad(1.0));
}