#include "gtest/gtest.h"

#include "PointSpreadFunction.h"

using namespace std;

TEST(PointSpreadFunctionTest, Constructor_ConfigurationParameters)
{
    ConfigurationParameters cp = ConfigurationParameters("../testData/input_PointSpreadFunctionTest.yaml");

    PointSpreadFunction psf = PointSpreadFunction(cp);

}