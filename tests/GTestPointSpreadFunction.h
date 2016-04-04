#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "PointSpreadFunction.h"
#include "Units.h"
#include "Constants.h"
#include "ArrayOperations.h"

using namespace std;


/**
 * \class PointSpreadFunctionTest
 * 
 * \brief Test Fixture for the PointSpreadFunction class. 
 * 
 * \details
 * 
 * Setup configuration parameters for the PointSpreadFunction class and handle the HDF5 open-close.
 * 
 * Input parameters for this test are located in the file 'input_PointSpreadFunctionTest.yaml' in the
 * testData directory of the distribution.
 * 
 * The test creates an HDF5 file 'pointSpreadFunctionTest.hdf5' in the current working directory. 
 * This file is removed after the test finishes.
 */

class PointSpreadFunctionTest : public testing::Test
{
    protected:

        virtual void SetUp()
        {
            configParams = ConfigurationParameters("../testData/input_PointSpreadFunctionTest.yaml");
        
            hdf5File.open(hdf5Filename);
        }

        virtual void TearDown()
        {
            hdf5File.close();
            FileUtilities::remove(hdf5Filename);
        }

        string hdf5Filename = "pointSpreadFunctionTest.hdf5";
        ConfigurationParameters configParams;
        HDF5File hdf5File;
};





TEST_F(PointSpreadFunctionTest, Constructor_ConfigurationParameters)
{

    LOG_STARTING_OF_TEST

    PointSpreadFunction psf = PointSpreadFunction(configParams, hdf5File);

}

TEST_F(PointSpreadFunctionTest, Selection)
{

    LOG_STARTING_OF_TEST

    PointSpreadFunction psf = PointSpreadFunction(configParams, hdf5File);
    psf.select(deg2rad(13.0));

}