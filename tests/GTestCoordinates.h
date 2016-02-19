#include "gtest/gtest.h"

#include "Units.h"
#include "Constants.h"
#include "Coordinates.h"

using namespace std;

// TODO: Find some good boundary tests for calculating the angular distance bewteen two coordinates


TEST(CoordinatesTest, angularDistance)
{
        Coordinates opticalAxis(0.0, 0.0, Angle::degrees);
        Coordinates star(178, -70, Angle::degrees);

        double angle = angularDistanceBetween(opticalAxis, star, Angle::degrees);
        Log.debug("CoordinatesTest.angularDistance: angle = " + to_string(angle));

}