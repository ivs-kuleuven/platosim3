#include <cstdio>
#include <cmath>
#include <map>

#include "gtest/gtest.h"

#include "Simulation.h"


using namespace std;

// Test Fixture: setup configuration parameters for the Detector class
//               and handle the HDF5 open-close.
//               
class DetectorTest : public testing::Test
{
    protected:

        virtual void SetUp()
        {
            cp_ = ConfigurationParameters("../testData/input_DetectorTest.yaml");

            remove(hdf5Filename.c_str());
            hdf5File_.open(hdf5Filename);
        }

        virtual void TearDown()
        {
            hdf5File_.close();
        }

        void reset()
        {
            hdf5File_.close();
            remove(hdf5Filename.c_str());
            hdf5File_.open(hdf5Filename);
        }

        string hdf5Filename = "/tmp/detectorTest.hdf5";
        ConfigurationParameters cp_;
        HDF5File hdf5File_;
};

// This subclass of Detector serves the sole purpose of testing protected methods of Detector.
// 
class MyDetector : public Detector
{
    public:
        MyDetector(ConfigurationParameters &configParam, HDF5File &hdf5File, Camera &camera);

        pair<double, double> test_pixelToFocalPlaneCoordinates(double row, double column) {return pixelToFocalPlaneCoordinates(row, column);};
        pair<double, double> test_focalPlaneToPixelCoordinates(double xFPprime, double yFPprime) {return focalPlaneToPixelCoordinates(xFPprime, yFPprime);};

};

MyDetector::MyDetector(ConfigurationParameters &configParam, HDF5File &hdf5File, Camera &camera)
: Detector(configParam, hdf5File, camera)
{
}


TEST_F(DetectorTest, checkConversionsBetweenPixelsAndFocalPlane)
{
    double row, column;
    double xFPprime, yFPprime;

    vector<map<string, double>> pixel2fp;

    JitterGenerator *jitterGenerator = new JitterFromRedNoise(cp_);
    Platform platform(cp_, hdf5File_, *jitterGenerator);
    Sky sky(cp_);
    Telescope telescope(cp_, hdf5File_, platform);
    Camera camera(cp_, hdf5File_, telescope, sky);

    // Settings for camera A
    pixel2fp.push_back(map<string, double> {{"xCCD",    0.0}, {"yCCD",    0.0}, {"xFP",  -1.0000}, {"yFP",  82.1620}, {"zeroPointX", -1.0000}, {"zeroPointY", 82.1620}, {"ccdAngle", 180.0}});
    pixel2fp.push_back(map<string, double> {{"xCCD", 4509.0}, {"yCCD",    0.0}, {"xFP", -82.1620}, {"yFP",  82.1620}, {"zeroPointX", -1.0000}, {"zeroPointY", 82.1620}, {"ccdAngle", 180.0}});
    pixel2fp.push_back(map<string, double> {{"xCCD", 4509.0}, {"yCCD", 4509.0}, {"xFP", -82.1620}, {"yFP",   1.0000}, {"zeroPointX", -1.0000}, {"zeroPointY", 82.1620}, {"ccdAngle", 180.0}});
    pixel2fp.push_back(map<string, double> {{"xCCD",    0.0}, {"yCCD", 4509.0}, {"xFP",  -1.0000}, {"yFP",   1.0000}, {"zeroPointX", -1.0000}, {"zeroPointY", 82.1620}, {"ccdAngle", 180.0}});

    // Settings for camera B
    pixel2fp.push_back(map<string, double> {{"xCCD",    0.0}, {"yCCD",    0.0}, {"xFP",  82.1620}, {"yFP",   1.0000}, {"zeroPointX", 82.1620}, {"zeroPointY",  1.0000}, {"ccdAngle",  90.0}});
    pixel2fp.push_back(map<string, double> {{"xCCD", 4509.0}, {"yCCD",    0.0}, {"xFP",  82.1620}, {"yFP",  82.1620}, {"zeroPointX", 82.1620}, {"zeroPointY",  1.0000}, {"ccdAngle",  90.0}});
    pixel2fp.push_back(map<string, double> {{"xCCD", 4509.0}, {"yCCD", 4509.0}, {"xFP",   1.0000}, {"yFP",  82.1620}, {"zeroPointX", 82.1620}, {"zeroPointY",  1.0000}, {"ccdAngle",  90.0}});
    pixel2fp.push_back(map<string, double> {{"xCCD",    0.0}, {"yCCD", 4509.0}, {"xFP",   1.0000}, {"yFP",   1.0000}, {"zeroPointX", 82.1620}, {"zeroPointY",  1.0000}, {"ccdAngle",  90.0}});

    for (auto &data: pixel2fp)
    {
        // FIXME: This is an inconvenience that the selection/orientation of the CCD can not be done with accessor methods.
        
        reset();
        cp_.setParameter("CCD/OriginOffsetX", to_string(data["zeroPointX"]));
        cp_.setParameter("CCD/OriginOffsetY", to_string(data["zeroPointY"]));
        cp_.setParameter("CCD/Orientation", to_string(data["ccdAngle"]));
    
        MyDetector detector = MyDetector(cp_, hdf5File_, camera);
    
        row = data["xCCD"];
        column = data["yCCD"];
        tie(xFPprime, yFPprime) = detector.test_pixelToFocalPlaneCoordinates(row, column);
    
        EXPECT_NEAR(data["xFP"], xFPprime, 0.00001);
        EXPECT_NEAR(data["yFP"], yFPprime, 0.00001);    
    }
}

