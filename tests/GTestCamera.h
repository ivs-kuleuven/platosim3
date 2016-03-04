#include <string>
#include <cstdio>
#include <map>

#include "gtest/gtest.h"
#include "gtest_definitions.h"

#include "Simulation.h"
#include "SkyCoordinates.h"
#include "StringUtilities.h"
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
        MyCamera(ConfigurationParameters &configParam, HDF5File &hdf5file, Telescope &telescope, Sky &sky): Camera(configParam, hdf5file, telescope, sky) {};

        pair<double, double> test_skyToAngularFocalPlaneCoordinates(double raStar, double decStar) {return skyToAngularFocalPlaneCoordinates(raStar, decStar);};
        pair<double, double> test_angularFocalPlaneToSkyCoordinates(double xFPprime, double yFPprime) {return angularFocalPlaneToSkyCoordinates(xFPprime, yFPprime);};

        pair<double, double> test_angularToPlanarFocalPlaneCoordinates(double xFPrad, double yFPrad) {return angularToPlanarFocalPlaneCoordinates(xFPrad, yFPrad);};
        pair<double, double> test_planarToAngularFocalPlaneCoordinates(double xFPmm, double yFPmm) {return planarToAngularFocalPlaneCoordinates(xFPmm, yFPmm);};

        pair<double, double> test_planarToDistortedFocalPlaneCoordinates(double xFPmm, double yFPmm) {return planarToDistortedFocalPlaneCoordinates(xFPmm, yFPmm);};

        double test_getGnomonicRadialDistanceFromOpticalAxis(double xFPprime, double yFPprime) {return getGnomonicRadialDistanceFromOpticalAxis(xFPprime, yFPprime);};

        void test_setDistortionPolynomial(Polynomial1D &polynomial) {setDistortionPolynomial(polynomial);};
};








// The table below represents the data from the list file as provided with the PSFs.
//
// The relation between xDeg, yDeg and radius is provided by the getGnomonicRadialDistanceFromOpticalAxisNormalized(..)
// method in Camera. The relation between xDeg, yDeg and xFP, yFP is currently not available and can not be tested.
// 
TEST_F(CameraTest, GnomonicRadialDistance)
{

    LOG_STARTING_OF_TEST

    using StringUtilities::dtos;

    // The values in this table come from different sources
    // xDeg, yDeg, and radius are taken from the table / attributes provided with the PSF delivery of <DATE>
    // xFP, yFP are the paraxial x and y field coordinates taken from the distorion table calculated from ZEMAX <DATE>

    vector<map<string, double>> gnomonic;

    gnomonic.push_back(map<string, double> {{"xDeg",  0.0}, {"yDeg",  0.0}, {"xFP",   0.000000}, {"yFP",   0.00000}, {"radius",   0.0000}});
    gnomonic.push_back(map<string, double> {{"xDeg",  1.0}, {"yDeg",  1.0}, {"xFP",   4.313600}, {"yFP",  4.313600}, {"radius",   1.4141}});
    gnomonic.push_back(map<string, double> {{"xDeg",  2.0}, {"yDeg",  2.0}, {"xFP",   8.629828}, {"yFP",  8.629828}, {"radius",   2.8273}});
    gnomonic.push_back(map<string, double> {{"xDeg",  3.0}, {"yDeg",  3.0}, {"xFP",  12.951322}, {"yFP", 12.951322}, {"radius",   4.2388}});
    gnomonic.push_back(map<string, double> {{"xDeg",  4.0}, {"yDeg",  4.0}, {"xFP",  17.280730}, {"yFP", 17.280730}, {"radius",   5.6477}});
    gnomonic.push_back(map<string, double> {{"xDeg",  5.0}, {"yDeg",  5.0}, {"xFP",  21.620719}, {"yFP", 21.620719}, {"radius",   7.0532}});
    gnomonic.push_back(map<string, double> {{"xDeg",  6.0}, {"yDeg",  6.0}, {"xFP",  25.973984}, {"yFP", 25.973984}, {"radius",   8.4545}});
    gnomonic.push_back(map<string, double> {{"xDeg",  7.0}, {"yDeg",  7.0}, {"xFP",  30.343251}, {"yFP", 30.343251}, {"radius",   9.8508}});
    gnomonic.push_back(map<string, double> {{"xDeg",  8.0}, {"yDeg",  8.0}, {"xFP",  34.731287}, {"yFP", 34.731287}, {"radius",  11.2413}});
    gnomonic.push_back(map<string, double> {{"xDeg",  9.0}, {"yDeg",  9.0}, {"xFP",  39.140906}, {"yFP", 39.140906}, {"radius",  12.6253}});
    gnomonic.push_back(map<string, double> {{"xDeg", 10.0}, {"yDeg", 10.0}, {"xFP",  43.574973}, {"yFP", 43.574973}, {"radius",  14.0019}});
    gnomonic.push_back(map<string, double> {{"xDeg", 11.0}, {"yDeg", 11.0}, {"xFP",  48.036419}, {"yFP", 48.036419}, {"radius",  15.3707}});
    gnomonic.push_back(map<string, double> {{"xDeg", 12.0}, {"yDeg", 12.0}, {"xFP",  52.528243}, {"yFP", 52.528243}, {"radius",  16.7308}});
    gnomonic.push_back(map<string, double> {{"xDeg", 13.0}, {"yDeg", 13.0}, {"xFP",  57.053521}, {"yFP", 57.053521}, {"radius",  18.0817}});
    gnomonic.push_back(map<string, double> {{"xDeg", 13.6}, {"yDeg", 13.6}, {"xFP",  59.786060}, {"yFP", 59.786060}, {"radius",  18.8876}});

    JitterGenerator *jitterGenerator = new JitterFromRedNoise(cp_);
    Platform platform = Platform(cp_, hdf5File_, *jitterGenerator);
    Telescope telescope = Telescope(cp_, hdf5File_, platform);
    Sky sky = Sky(cp_);

    MyCamera camera = MyCamera(cp_, hdf5File_, telescope, sky);

    double raStar, decStar;      // [rad]
    double xFPprime, yFPprime;   // [mm]
    double xDeg, yDeg;           // [deg]
    double xFPrad, yFPrad;       // [rad]
    double radius;               // [rad]

    for (auto &data: gnomonic)
    {
        // Other tests that need to be performed is the conversion from xFP, yFP to xDeg, yDeg (in the data table above)

        xDeg = data["xDeg"];
        yDeg = data["yDeg"];
        xFPrad = deg2rad(xDeg);
        yFPrad = deg2rad(yDeg);

        radius = camera.test_getGnomonicRadialDistanceFromOpticalAxis(xFPrad, yFPrad); // [radians] -> [radians]

        EXPECT_NEAR(data["radius"], rad2deg(radius), 0.0001);

        Log.debug("CameraTest.GnomonicRadialDistance: xDeg, yDeg [deg] = " + dtos(xDeg) + ", " + dtos(yDeg));
        Log.debug("CameraTest.GnomonicRadialDistance: xFPrad, yFPrad [rad] = " + dtos(xFPrad) + ", " + dtos(yFPrad));
        Log.debug("CameraTest.GnomonicRadialDistance: radius [rad] = " + dtos(radius));
        Log.debug("CameraTest.GnomonicRadialDistance: radius [deg] = " + dtos(rad2deg(radius)));

        tie(xFPprime, yFPprime) = camera.test_angularToPlanarFocalPlaneCoordinates(xFPrad, yFPrad);

        EXPECT_NEAR(data["xFP"], xFPprime, 0.00001);
        EXPECT_NEAR(data["yFP"], yFPprime, 0.00001);

        Log.debug("CameraTest.GnomonicRadialDistance: xFPprime, yFPprime [mm]= " + dtos(xFPprime) + ", " + dtos(yFPprime));

        tie(raStar, decStar) = camera.test_angularFocalPlaneToSkyCoordinates(xFPrad, yFPrad);

        Log.debug("CameraTest.GnomonicRadialDistance: raStar, decStar [rad] = " + dtos(raStar) + ", " + dtos(decStar));
        Log.debug("CameraTest.GnomonicRadialDistance: raStar, decStar [deg] = " + dtos(rad2deg(raStar)) + ", " + dtos(rad2deg(decStar)));

    }

}



TEST_F(CameraTest, distortedCoordinates)
{
    LOG_STARTING_OF_TEST

    using StringUtilities::dtos;

    JitterGenerator *jitterGenerator = new JitterFromRedNoise(cp_);
    Platform platform = Platform(cp_, hdf5File_, *jitterGenerator);
    Telescope telescope = Telescope(cp_, hdf5File_, platform);
    Sky sky = Sky(cp_);

    int degree = 1;
    vector<double> coeff = {0.0, 1.0};
    Polynomial1D polynomial = Polynomial1D(degree, coeff);

    MyCamera camera = MyCamera(cp_, hdf5File_, telescope, sky);

    camera.test_setDistortionPolynomial(polynomial);

    double xFPdist, yFPdist;   // [mm]

    tie(xFPdist, yFPdist) = camera.test_planarToDistortedFocalPlaneCoordinates(10.0, 0.0);
    EXPECT_NEAR(10.0000, xFPdist, 0.00001);
    EXPECT_NEAR( 0.0000, yFPdist, 0.00001);

    tie(xFPdist, yFPdist) = camera.test_planarToDistortedFocalPlaneCoordinates(0.0, 10.0);
    EXPECT_NEAR( 0.0000, xFPdist, 0.00001);
    EXPECT_NEAR(10.0000, yFPdist, 0.00001);

    tie(xFPdist, yFPdist) = camera.test_planarToDistortedFocalPlaneCoordinates(5.0, 5.0);
    EXPECT_NEAR( 5.0000, xFPdist, 0.00001);
    EXPECT_NEAR( 5.0000, yFPdist, 0.00001);

    degree = 2;
    coeff = {2.0, 0.5, 1.5};
    polynomial = Polynomial1D(degree, coeff);

    camera.test_setDistortionPolynomial(polynomial);

    tie(xFPdist, yFPdist) = camera.test_planarToDistortedFocalPlaneCoordinates(10.0, 0.0);
    EXPECT_NEAR(157.0000, xFPdist, 0.00001);
    EXPECT_NEAR( 2.0000, yFPdist, 0.00001);

    tie(xFPdist, yFPdist) = camera.test_planarToDistortedFocalPlaneCoordinates(0.0, 10.0);
    EXPECT_NEAR( 2.0000, xFPdist, 0.00001);
    EXPECT_NEAR(157.0000, yFPdist, 0.00001);

    tie(xFPdist, yFPdist) = camera.test_planarToDistortedFocalPlaneCoordinates(5.0, 5.0);
    EXPECT_NEAR( 42.0000, xFPdist, 0.00001);
    EXPECT_NEAR( 42.0000, yFPdist, 0.00001);

    Log.debug("CameraTest.distortedCoordinates: xFPdist, yFPdist = " + dtos(xFPdist) + ", " + dtos(yFPdist));


}


TEST_F(CameraTest, reproduceDistortionMap)
{
    LOG_STARTING_OF_TEST

    using StringUtilities::dtos;

    JitterGenerator *jitterGenerator = new JitterFromRedNoise(cp_);
    Platform platform = Platform(cp_, hdf5File_, *jitterGenerator);
    Telescope telescope = Telescope(cp_, hdf5File_, platform);
    Sky sky = Sky(cp_);

    // Just a few values from the current distortion map for which the Polynomial1D was fitted.
    vector<map<string, double>> distortion;
    distortion.push_back(map<string, double> {{"xFPmm",  0.000000}, {"xFPdist",   0.00000}});
    distortion.push_back(map<string, double> {{"xFPmm", 10.789752}, {"xFPdist",  10.796226}});
    distortion.push_back(map<string, double> {{"xFPmm", 24.666449}, {"xFPdist",  24.743989}});
    distortion.push_back(map<string, double> {{"xFPmm", 68.998805}, {"xFPdist",  70.734282}});
    distortion.push_back(map<string, double> {{"xFPmm", 80.296089}, {"xFPdist",  83.062336}});

    // These values are for a fit of Polynomial1D with degree=3 to the distortion table

    int degree = 3;
    vector<double> coeff = {-0.0036696919678, 1.0008542317, -4.12553764817e-05, 5.7201219949e-06};
    Polynomial1D polynomial = Polynomial1D(degree, coeff);

    MyCamera camera = MyCamera(cp_, hdf5File_, telescope, sky);
    camera.test_setDistortionPolynomial(polynomial);

    double xFPdist, yFPdist;   // [mm]

    for (auto &data: distortion)
    {
        double xFPmm = data["xFPmm"];

        tie(xFPdist, yFPdist) = camera.test_planarToDistortedFocalPlaneCoordinates(xFPmm, xFPmm);
        EXPECT_NEAR( data["xFPdist"], xFPdist, 0.006);
        EXPECT_NEAR( data["xFPdist"], yFPdist, 0.006);
    }

    // These values are for a fit of Polynomial1D with degree=4 to the distortion table
    // Obviously the error of the fit is smaller, but we need to check what the effect is on the noise level.
    // Remember that 1 pixel = 0.018 mm

    degree = 4;
    coeff = {0.000814670391532, 0.99970324817, 2.39592182367e-05, 4.44973838376e-06, 7.93413878401e-09};
    polynomial = Polynomial1D(degree, coeff);

    camera.test_setDistortionPolynomial(polynomial);

    for (auto &data: distortion)
    {
        double xFPmm = data["xFPmm"];

        tie(xFPdist, yFPdist) = camera.test_planarToDistortedFocalPlaneCoordinates(xFPmm, xFPmm);
        EXPECT_NEAR( data["xFPdist"], xFPdist, 0.002);
        EXPECT_NEAR( data["xFPdist"], yFPdist, 0.002);
    }

}




