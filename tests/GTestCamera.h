#include <string>
#include <cstdio>
#include <map>

#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "Simulation.h"
#include "SkyCoordinates.h"
#include "FileUtilities.h"

/**
 * \class CameraTest
 * 
 * \brief Test Fixture for the Camera class. 
 * 
 * \details
 * 
 * Setup configuration parameters for the Camera class and handle the HDF5 open-close.
 * 
 * Input parameters for this test are located in the file 'input_CameraTest.yaml' in the
 * testData directory of the distribution.
 * 
 * The test creates an HDF5 file 'cameraTest.hdf5' in the current working directory. This file is removed after the test finishes.
 */

class CameraTest : public testing::Test
{
    protected:

        virtual void SetUp()
        {
            cp_ = ConfigurationParameters("../testData/input_CameraTest.yaml");
        
            hdf5File_.open(hdf5Filename);
        }

        virtual void TearDown()
        {
            hdf5File_.close();
            FileUtilities::remove(hdf5Filename);
        }

        string hdf5Filename = "cameraTest.hdf5";
        ConfigurationParameters cp_;
        HDF5File hdf5File_;
};




/**
 * 
 * \brief This subclass of Camera serves the sole purpose of testing protected methods of Camera.
 * 
 */

class MyCamera : public Camera
{
    public:
        MyCamera(ConfigurationParameters &configParam, HDF5File &hdf5file, Telescope &telescope, Sky &sky);

        pair<double, double> test_skyToFocalPlaneCoordinates(double raStar, double decStar) {return skyToFocalPlaneCoordinates(raStar, decStar);};
        pair<double, double> test_focalPlaneToSkyCoordinates(double x, double y) {return focalPlaneToSkyCoordinates(x, y);};
        double test_getGnomonicRadialDistanceFromOpticalAxis(double xFPprime, double yFPprime) {return getGnomonicRadialDistanceFromOpticalAxis(xFPprime, yFPprime);};
        double test_getGnomonicRadialDistanceFromOpticalAxisNormalized(double xFPprime, double yFPprime) {return getGnomonicRadialDistanceFromOpticalAxisNormalized(xFPprime, yFPprime);};
};


MyCamera::MyCamera(ConfigurationParameters &configParam, HDF5File &hdf5file, Telescope &telescope, Sky &sky)
: Camera(configParam, hdf5file, telescope, sky)
{
}







// The table below represents the data from the list file as provided with the PSFs.
//
// The relation between xDeg, yDeg and radius is provided by the getGnomonicRadialDistanceFromOpticalAxisNormalized(..)
// method in Camera. The relation between xDeg, yDeg and xFP, yFP is currently not available and can not be tested.
// 
TEST_F(CameraTest, GnomonicRadialDistance) {

    LOG_STARTING_OF_TEST

    vector<map<string, double>> gnomonic;

    gnomonic.push_back(map<string, double> {{"xDeg",  0.0}, {"yDeg",  0.0}, {"xFP",   0.0000}, {"yFP",   0.0000}, {"radius",   0.0000}});
    gnomonic.push_back(map<string, double> {{"xDeg",  1.0}, {"yDeg",  1.0}, {"xFP",   0.0000}, {"yFP",   0.0000}, {"radius",   1.4141}});
    gnomonic.push_back(map<string, double> {{"xDeg",  2.0}, {"yDeg",  2.0}, {"xFP",   8.6354}, {"yFP",   8.6354}, {"radius",   2.8273}});
    gnomonic.push_back(map<string, double> {{"xDeg",  3.0}, {"yDeg",  3.0}, {"xFP",   0.0000}, {"yFP",   0.0000}, {"radius",   4.2388}});
    gnomonic.push_back(map<string, double> {{"xDeg",  4.0}, {"yDeg",  4.0}, {"xFP",  17.3325}, {"yFP",  17.3325}, {"radius",   5.6477}});
    gnomonic.push_back(map<string, double> {{"xDeg",  5.0}, {"yDeg",  5.0}, {"xFP",   0.0000}, {"yFP",   0.0000}, {"radius",   7.0532}});
    gnomonic.push_back(map<string, double> {{"xDeg",  6.0}, {"yDeg",  6.0}, {"xFP",  26.1552}, {"yFP",  26.1552}, {"radius",   8.4545}});
    gnomonic.push_back(map<string, double> {{"xDeg",  7.0}, {"yDeg",  7.0}, {"xFP",  30.6341}, {"yFP",  30.6341}, {"radius",   9.8508}});
    gnomonic.push_back(map<string, double> {{"xDeg",  8.0}, {"yDeg",  8.0}, {"xFP",  35.1698}, {"yFP",  35.1698}, {"radius",  11.2413}});
    gnomonic.push_back(map<string, double> {{"xDeg",  9.0}, {"yDeg",  9.0}, {"xFP",   0.0000}, {"yFP",   0.0000}, {"radius",  12.6253}});
    gnomonic.push_back(map<string, double> {{"xDeg", 10.0}, {"yDeg", 10.0}, {"xFP",  44.4487}, {"yFP",  44.4487}, {"radius",  14.0019}});
    gnomonic.push_back(map<string, double> {{"xDeg", 11.0}, {"yDeg", 11.0}, {"xFP",   0.0000}, {"yFP",   0.0000}, {"radius",  15.3707}});
    gnomonic.push_back(map<string, double> {{"xDeg", 12.0}, {"yDeg", 12.0}, {"xFP",  54.0726}, {"yFP",  54.0726}, {"radius",  16.7308}});
    gnomonic.push_back(map<string, double> {{"xDeg", 13.0}, {"yDeg", 13.0}, {"xFP",   0.0000}, {"yFP",   0.0000}, {"radius",  18.0817}});
    gnomonic.push_back(map<string, double> {{"xDeg", 13.6}, {"yDeg", 13.6}, {"xFP",  62.0957}, {"yFP",  62.0957}, {"radius",  18.8876}});

    JitterGenerator *jitterGenerator = new JitterFromRedNoise(cp_);
    Platform platform = Platform(cp_, hdf5File_, *jitterGenerator);
    Telescope telescope = Telescope(cp_, hdf5File_, platform);
    Sky sky = Sky(cp_);

    MyCamera camera = MyCamera(cp_, hdf5File_, telescope, sky);

    double xFPprime, yFPprime;   // [rad]
    double radius;               // [rad]

    for (auto &data: gnomonic)
    {
        // Other tests that need to be performed is the conversion from xFP, yFP to xDeg, yDeg (in the data table above)
        
        // xFPprime and yFPprime as described in eq. 24 of TN on PlatSim Reference Frames v1.3
        xFPprime = deg2rad(data["xDeg"]);
        yFPprime = deg2rad(data["yDeg"]);

        radius = camera.test_getGnomonicRadialDistanceFromOpticalAxisNormalized(xFPprime, yFPprime); // [radians] -> [radians]

        EXPECT_NEAR(data["radius"], rad2deg(radius), 0.0001);

        Log.debug("CameraTest.GnomonicRadialDistance: xFPprime, yFPprime = " + to_string(xFPprime) + ", " + to_string(yFPprime));
        Log.debug("CameraTest.GnomonicRadialDistance: radius = " + to_string(radius));

    }

}

