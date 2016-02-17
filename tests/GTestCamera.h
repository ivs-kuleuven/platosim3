#include <cstdio>
#include <map>

#include "gtest/gtest.h"

#include "Simulation.h"
#include "Coordinates.h"


// Test Fixture: setup configuration parameters for the Camera class
//               and handle the HDF5 open-close.
//               
class CameraTest : public testing::Test
{
    protected:

        virtual void SetUp()
        {
            cp_ = ConfigurationParameters("../testData/input_CameraTest.yaml");
        
            remove(hdf5Filename.c_str());
            hdf5File_.open(hdf5Filename);
        }

        virtual void TearDown()
        {
            hdf5File_.close();
        }

        string hdf5Filename = "/tmp/cameraTest.hdf5";
        ConfigurationParameters cp_;
        HDF5File hdf5File_;
};

// This subclass of Camera serves the sole purpose of testing protected methods of Camera.
// 
class MyCamera : public Camera
{
    public:
        MyCamera(ConfigurationParameters &configParam, HDF5File &hdf5file, Telescope &telescope, Sky &sky);

        pair<double, double> test_skyToFocalPlaneCoordinates(double raStar, double decStar) {return skyToFocalPlaneCoordinates(raStar, decStar);};
        pair<double, double> test_focalPlaneToSkyCoordinates(double x, double y) {return focalPlaneToSkyCoordinates(x, y);};
        double test_getGnomonicRadialDistance(double xDeg, double yDeg) {return getGnomonicRadialDistanceFromOpticalAxis(xDeg, yDeg);};
};


MyCamera::MyCamera(ConfigurationParameters &configParam, HDF5File &hdf5file, Telescope &telescope, Sky &sky)
: Camera(configParam, hdf5file, telescope, sky)
{
}


// This test is DISABLED_ in the public master as I'm still working on the code and the tests fail currently.
// The code itself is not activated yet in platosim.
TEST_F(CameraTest, DISABLED_skyToFocalPlaneCoordinates) {

    vector<map<string, double>> gnomonic;

    gnomonic.push_back(map<string, double> {{"xDeg",  0.0}, {"yDeg",  0.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",   0.0000}});
    gnomonic.push_back(map<string, double> {{"xDeg",  1.0}, {"yDeg",  1.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",   1.4141}});
    gnomonic.push_back(map<string, double> {{"xDeg",  2.0}, {"yDeg",  2.0}, {"xFP",  8.6354}, {"yFP",  8.6354}, {"radialDistance",   2.8273}});
    gnomonic.push_back(map<string, double> {{"xDeg",  3.0}, {"yDeg",  3.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",   4.2388}});
    gnomonic.push_back(map<string, double> {{"xDeg",  4.0}, {"yDeg",  4.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",   5.6477}});
    gnomonic.push_back(map<string, double> {{"xDeg",  5.0}, {"yDeg",  5.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",   7.0532}});
    gnomonic.push_back(map<string, double> {{"xDeg",  6.0}, {"yDeg",  6.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",   8.4545}});
    gnomonic.push_back(map<string, double> {{"xDeg",  7.0}, {"yDeg",  7.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",   9.8508}});
    gnomonic.push_back(map<string, double> {{"xDeg",  8.0}, {"yDeg",  8.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",  11.2413}});
    gnomonic.push_back(map<string, double> {{"xDeg",  9.0}, {"yDeg",  9.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",  12.6253}});
    gnomonic.push_back(map<string, double> {{"xDeg", 10.0}, {"yDeg", 10.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",  14.0019}});
    gnomonic.push_back(map<string, double> {{"xDeg", 11.0}, {"yDeg", 11.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",  15.3707}});
    gnomonic.push_back(map<string, double> {{"xDeg", 12.0}, {"yDeg", 12.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",  16.7308}});
    gnomonic.push_back(map<string, double> {{"xDeg", 13.0}, {"yDeg", 13.0}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",  18.0817}});
    gnomonic.push_back(map<string, double> {{"xDeg", 13.6}, {"yDeg", 13.6}, {"xFP",  0.0000}, {"yFP",  0.0000}, {"radialDistance",  18.8876}});

    JitterGenerator *jitterGenerator = new JitterFromRedNoise(cp_);
    Platform platform = Platform(cp_, hdf5File_, *jitterGenerator);
    Telescope telescope = Telescope(cp_, hdf5File_, platform);
    Sky sky = Sky(cp_);

    MyCamera camera = MyCamera(cp_, hdf5File_, telescope, sky);

    // TODO:
    // Test incomplete, I now only test if the from and to methods end up with the original values.
    // What should be done is check for boundary values if the expected result is returned.
    // Other tests should be done with optical axis not being (0,0).

    double xFPprime, yFPprime;   // [mm]
    double raStar, decStar;      // [rad]
    double raStar2, decStar2;    // [rad]
    double radialDistance;       // [deg]

    for (auto &data: gnomonic)
    {
        raStar = deg2rad(data["raStar"]);
        decStar = deg2rad(data["decStar"]);

        tie(xFPprime, yFPprime) = camera.test_skyToFocalPlaneCoordinates(raStar, decStar);    // [radians] -> [mm]

        tie(raStar2, decStar2) = camera.test_focalPlaneToSkyCoordinates(xFPprime, yFPprime);  // [mm] -> [radians]

        radialDistance = camera.test_getGnomonicRadialDistance(xFPprime, yFPprime); // [mm] -> [degrees]

        EXPECT_DOUBLE_EQ(raStar, raStar2);
        EXPECT_DOUBLE_EQ(decStar, decStar2);

        //EXPECT_NEAR(data["xFPprime"], xFPprime, 0.0001);
        //EXPECT_NEAR(data["yFPprime"], yFPprime, 0.0001);
        EXPECT_NEAR(data["radialDistance"], radialDistance, 0.00001);


        Log.debug("CameraTest.skyToFocalPlaneCoordinates: raStar, decStar = " + to_string(raStar) + ", " + to_string(decStar));
        Log.debug("CameraTest.skyToFocalPlaneCoordinates: raStar2, decStar2 = " + to_string(raStar2) + ", " + to_string(decStar2));
        Log.debug("CameraTest.skyToFocalPlaneCoordinates: xFPprime, yFPprime = " + to_string(xFPprime) + ", " + to_string(yFPprime));
        Log.debug("CameraTest.skyToFocalPlaneCoordinates: radialDistance = " + to_string(radialDistance));

    }

    tie(raStar, decStar) = camera.test_focalPlaneToSkyCoordinates(8.6354, 8.6354);

    Coordinates opticalAxis(0.0, 0.0, Angle::degrees);
    Coordinates starPosition(raStar, decStar, Angle::radians);

    double ad = angularDistanceBetween(opticalAxis, starPosition, Angle::degrees);
    Log.debug("CameraTest.skyToFocalPlaneCoordinates: 8.6354, 8.6354 -> angularDistance = " + to_string(ad));
    Log.debug("CameraTest.skyToFocalPlaneCoordinates: 8.6354, 8.6354 -> raStar, decStar = " + to_string(rad2deg(raStar)) + ", " + to_string(rad2deg(decStar)));

    tie(raStar, decStar) = camera.test_focalPlaneToSkyCoordinates(17.3325, 17.3325);

    Coordinates starPosition2(raStar, decStar, Angle::radians);

    double ad2 = angularDistanceBetween(opticalAxis, starPosition2, Angle::degrees);
    Log.debug("CameraTest.skyToFocalPlaneCoordinates: 17.3325, 17.3325 -> angularDistance = " + to_string(ad2));
    Log.debug("CameraTest.skyToFocalPlaneCoordinates: 17.3325, 17.3325 -> raStar, decStar = " + to_string(rad2deg(raStar)) + ", " + to_string(rad2deg(decStar)));

}

