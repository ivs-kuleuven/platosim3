#include "gtest/gtest.h"

#include "PointSpreadFunction.h"
#include "Units.h"
#include "Constants.h"
#include "ArrayOperations.h"

using namespace std;


TEST(PointSpreadFunctionTest, Constructor_ConfigurationParameters)
{
    ConfigurationParameters cp = ConfigurationParameters("../testData/input_PointSpreadFunctionTest.yaml");

    PointSpreadFunction psf = PointSpreadFunction(cp);

}

