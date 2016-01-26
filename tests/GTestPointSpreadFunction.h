#include "gtest/gtest.h"

#include "PointSpreadFunction.h"

using namespace std;

TEST(PointSpreadFunctionTest, Constructor_ConfigurationParameters)
{
    ConfigurationParameters cp = ConfigurationParameters("../testData/input.yaml");

    // This makes sure the HDF5 file is loaded from the correct location.
    cp.setNode("General/ProjectLocation", "../");

    PointSpreadFunction psf = PointSpreadFunction(cp);

}