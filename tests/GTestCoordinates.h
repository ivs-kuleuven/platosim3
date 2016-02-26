#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "Units.h"
#include "Constants.h"
#include "SkyCoordinates.h"

using namespace std;

// TODO: Find some good boundary tests for calculating the angular distance bewteen two coordinates


TEST(CoordinatesTest, angularDistance)
{
    LOG_STARTING_OF_TEST

    SkyCoordinates opticalAxis(0.0, 0.0, Angle::degrees);
    SkyCoordinates star(178, -70, Angle::degrees);

    double angle = angularDistanceBetween(opticalAxis, star, Angle::degrees);
    Log.debug("CoordinatesTest.angularDistance: angle = " + to_string(angle));

}