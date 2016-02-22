#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "PointSpreadFunction.h"
#include "Units.h"
#include "Constants.h"
#include "ArrayOperations.h"

using namespace std;


TEST(PointSpreadFunctionTest, Constructor_ConfigurationParameters)
{

    LOG_STARTING_OF_TEST

    ConfigurationParameters cp = ConfigurationParameters("../testData/input_PointSpreadFunctionTest.yaml");

    PointSpreadFunction psf = PointSpreadFunction(cp);

}

TEST(PointSpreadFunctionTest, Selection)
{

    LOG_STARTING_OF_TEST

    ConfigurationParameters cp = ConfigurationParameters("../testData/input_PointSpreadFunctionTest.yaml");

    PointSpreadFunction psf = PointSpreadFunction(cp);
    psf.select(deg2rad(13.0));

}