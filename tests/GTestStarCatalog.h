#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "Units.h"
#include "Constants.h"
#include "StarCatalog.h"
#include "Simulation.h"

using namespace std;


/**
 * \class StarCatalogTest
 * 
 * \brief Test Fixture for the StarCatalog class. 
 * 
 * \details
 * 
 * Setup configuration parameters for the StarCatalog class and handle the HDF5 open-close.
 * 
 * Input parameters for this test are located in the file 'input_StarCatalogTest.yaml' in the
 * testData directory of the distribution.
 * 
 * The test creates an HDF5 file 'StarCatalogTest.hdf5' in the current working directory. 
 * This file is removed after the test finishes.
 */

class StarCatalogTest : public testing::Test
{
    protected:

        virtual void SetUp()
        {
            configParams = ConfigurationParameters("../testData/input_StarCatalogTest.yaml");
        
            hdf5File.open(hdf5Filename);
        }

        virtual void TearDown()
        {
            hdf5File.close();
            FileUtilities::remove(hdf5Filename);
        }

        string hdf5Filename = "starCatalogTest.hdf5";
        ConfigurationParameters configParams;
        HDF5File hdf5File;
};


// Test the aberration correction for a one year observation


TEST_F(StarCatalogTest, aberrate)
{
    LOG_STARTING_OF_TEST

    StarCatalog starCatalog;

    // Add a few stars to the star catalog.

    starCatalog.addStar(0, 178.8, -70.9, 14.099, Angle::degrees);
    starCatalog.addStar(1, 163.8,  65.9, 14.099, Angle::degrees);
    starCatalog.addStar(2, 45.0,   45.0, 14.099, Angle::degrees);

    EXPECT_EQ(3, starCatalog.size());

    // Use the following variables for calculating the timeMiddle.
    // The observation takes one year, with 1 day exposures.

    int numExposures =  365;
    int exposureTime =  60*60*24;
    int readoutTime  =  0;

    // starIndex points into the starCatalog to choose a particular star

    int starIndex    =  2;

    // Set a proper platform pointing with respect to the used star coordinates

    switch(starIndex)
    {
        case 0:
            configParams.setParameter("ObservingParameters/RApointing",  "180");
            configParams.setParameter("ObservingParameters/DecPointing", "-70");
            break;
        case 1:
            configParams.setParameter("ObservingParameters/RApointing",  "163.9");
            configParams.setParameter("ObservingParameters/DecPointing", "65.9");
            break;
        case 2:
            configParams.setParameter("ObservingParameters/RApointing",  "45.0");
            configParams.setParameter("ObservingParameters/DecPointing", "54.0");
            break;
        default:
            FAIL() << "Incorrect startIndex given: " + to_string(starIndex);
    }

    // Remember the original coordinates of the star

    StarRecord starRecord = starCatalog[starIndex];
    double startRA  = starRecord.RA;
    double startDec = starRecord.dec;

    // Initialise the Platform (needed for raSun and decSun)

    JitterFromRedNoise jitterGenerator(configParams);    
    Platform platform = Platform(configParams, hdf5File, jitterGenerator);

    // The middle of the observation is defined as the point where the sun is oposite to the platform coordinates

    double timeMiddle = numExposures * (exposureTime + readoutTime) / 2.0;

    // Calculate and log the absolute aberrated coordinates for each day in one year

    string type = "absolute";

    for (double startTime = 0.0; startTime < timeMiddle * 2.0; startTime += exposureTime)
    {
        StarCatalog newStarCatalog = starCatalog.aberrate(platform, type, startTime, timeMiddle);
        StarRecord starRecord = newStarCatalog[starIndex];
        Log.debug("StarCatalogTest::aberrate: [" + type + "] " + to_string(rad2deg(starRecord.RA-startRA) * 3600.0 * cos(starRecord.dec)) + ", " + to_string(rad2deg(starRecord.dec-startDec) * 3600.0));
    }

    // Calculate and log the differential aberrated coordinates for each day in one year

    type = "differential";

    for (double startTime = 0.0; startTime < timeMiddle * 2.0; startTime += exposureTime)
    {
        StarCatalog newStarCatalog = starCatalog.aberrate(platform, "differential", startTime, timeMiddle);
        StarRecord starRecord = newStarCatalog[starIndex];
        Log.debug("StarCatalogTest::aberrate: [" + type + "] " + to_string(rad2deg(starRecord.RA-startRA) * 3600.0 * cos(starRecord.dec)) + ", " + to_string(rad2deg(starRecord.dec-startDec) * 3600.0));
    }

    // When you want to plot the absolute and differential aberration, this can be easily done in Python.
    // Therefore filter out the lines for the two types of aberration:
    //
    // $ grep "\[absolute\]" log.txt |cut -c 66- > absAberration.txt
    // $ grep "\[differential\]" log.txt |cut -c 70- > diffAberation.txt
    //
    // The use the code below to read the data and plot.
    //
    // import bokeh
    // from bokeh.plotting import figure, show
    // 
    // filenameDiff = "/Users/rik/Git/PlatoSim3/build/diffAberration.txt"
    // filenameAbs = "/Users/rik/Git/PlatoSim3/build/absAberration.txt"
    // 
    // with open(filenameDiff, 'r') as orbitFileDiff:
    //     raDiff, decDiff = [], []
    //     for line in orbitFileDiff.readlines():
    //         x, y = line.split(',')
    //         raDiff.append(float(x))
    //         decDiff.append(float(y))
    // 
    // 
    // with open(filenameAbs, 'r') as orbitFileAbs:
    //     raAbs, decAbs = [], []
    //     for line in orbitFileAbs.readlines():
    //         x, y = line.split(',')
    //         raAbs.append(float(x))
    //         decAbs.append(float(y))
    // 
    // 
    // fig = figure(width=1000, height=500, x_axis_label="(ra-ra0)*cos(dec) [arcsec]", y_axis_label="dec-dec0 [arcsec]")
    // fig.circle(raDiff, decDiff, size=2, color="red", alpha=0.5, legend="Differential")
    // fig.circle(raAbs, decAbs, size=2, color="navy", legend="Absolute")
    // show(fig)

}



