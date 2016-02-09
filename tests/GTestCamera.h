#include <map>

#include "gtest/gtest.h"


#define UNIT_TESTS

#include "Camera.h"
#include "Units.h"
#include "HDF5File.h"


TEST(CameraTest, skyToFocalPlaneCoordinates) {

    vector<map<string, double>> gnomonic;
    gnomonic.push_back(map<string, double> {{"raStar",  0.0}, {"decStar",  0.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",   0.0000}});
    gnomonic.push_back(map<string, double> {{"raStar",  1.0}, {"decStar",  1.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",   1.4141}});
    gnomonic.push_back(map<string, double> {{"raStar",  2.0}, {"decStar",  2.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",   2.8273}});
    gnomonic.push_back(map<string, double> {{"raStar",  3.0}, {"decStar",  3.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",   4.2388}});
    gnomonic.push_back(map<string, double> {{"raStar",  4.0}, {"decStar",  4.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",   5.6477}});
    gnomonic.push_back(map<string, double> {{"raStar",  5.0}, {"decStar",  5.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",   7.0532}});
    gnomonic.push_back(map<string, double> {{"raStar",  6.0}, {"decStar",  6.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",   8.4545}});
    gnomonic.push_back(map<string, double> {{"raStar",  7.0}, {"decStar",  7.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",   9.8508}});
    gnomonic.push_back(map<string, double> {{"raStar",  8.0}, {"decStar",  8.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",  11.2413}});
    gnomonic.push_back(map<string, double> {{"raStar",  9.0}, {"decStar",  9.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",  12.6253}});
    gnomonic.push_back(map<string, double> {{"raStar", 10.0}, {"decStar", 10.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",  14.0019}});
    gnomonic.push_back(map<string, double> {{"raStar", 11.0}, {"decStar", 11.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",  15.3707}});
    gnomonic.push_back(map<string, double> {{"raStar", 12.0}, {"decStar", 12.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",  16.7308}});
    gnomonic.push_back(map<string, double> {{"raStar", 13.0}, {"decStar", 13.0}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",  18.0817}});
    gnomonic.push_back(map<string, double> {{"raStar", 13.6}, {"decStar", 13.6}, {"x", 0.0}, {"y", 0.0}, {"radialCoordinate",  18.8876}});

    ConfigurationParameters cp = ConfigurationParameters();
    cp.setParameter("Camera/PlateScale", "0.833333333");
    cp.setParameter("Camera/FocalPlaneOrientation", "0.0");
    cp.setParameter("Camera/ThroughputBandwidth", "650");
    cp.setParameter("Camera/ThroughputLambdaC", "550");     
    cp.setParameter("ObservingParameters/RApointing", "0.0");
    cp.setParameter("ObservingParameters/DecPointing", "0.0");
    cp.setParameter("ObservingParameters/StarCatalogFile", "inputfiles/starcatalog.txt"); 
    cp.setParameter("Telescope/LightCollectingArea", "113.1");
    cp.setParameter("Telescope/TransmissionEfficiency", "0.757");
    cp.setParameter("Telescope/FOVSquareDegrees", "1072.0");  
    cp.setParameter("Telescope/DriftYawRms", "14.0");
    cp.setParameter("Telescope/DriftPitchRms", "2.3");       
    cp.setParameter("Telescope/DriftRollRms", "2.3");        
    cp.setParameter("Telescope/DriftTimeScale", "3600.");  

    HDF5File hdf5File;
    hdf5File.open("/tmp/cameraTest.hdf5");
    
    Telescope telescope = Telescope(cp, hdf5File);
    Sky sky = Sky(cp);

    Camera camera = Camera(cp, hdf5File, telescope, sky);

    // TODO:
    // Test incomplete, I now only test if the from and to methods end up with the original values.
    // What should be done is check for boundary values if the expected result is returned.
    // Other tests should be done with optical axis not being (0,0).

    double x, y;
    double raStar, decStar;
    double raStar2, decStar2;

    for (auto &data: gnomonic)
    {
        raStar = deg2rad(data["raStar"]);
        decStar = deg2rad(data["decStar"]);
        Log.debug("CameraTest.skyToFocalPlaneCoordinates: raStar, decStar = " + to_string(raStar) + ", " + to_string(decStar));

        tie(x, y) = camera.test_skyToFocalPlaneCoordinates(raStar, decStar);
        Log.debug("CameraTest.skyToFocalPlaneCoordinates: x, y = " + to_string(x) + ", " + to_string(y));

        tie(raStar2, decStar2) = camera.test_focalPlaneToSkyCoordinates(x, y);
        Log.debug("CameraTest.skyToFocalPlaneCoordinates: raStar2, decStar2 = " + to_string(raStar2) + ", " + to_string(decStar2));

        EXPECT_DOUBLE_EQ(raStar, raStar2);
        EXPECT_DOUBLE_EQ(decStar, decStar2);
    }

    hdf5File.close();
}

