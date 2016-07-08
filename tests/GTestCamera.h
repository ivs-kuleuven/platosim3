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
            configParams = ConfigurationParameters("../testData/input_CameraTest.yaml");
        
            hdf5File.open(hdf5Filename);
        }

        virtual void TearDown()
        {
            hdf5File.close();
            FileUtilities::remove(hdf5Filename);
        }

        string hdf5Filename = "cameraTest.hdf5";
        ConfigurationParameters configParams;
        HDF5File hdf5File;
};




/**
 * 
 * \brief This subclass of Camera serves the sole purpose of testing protected methods of Camera.
 * 
 */

class MyCamera : public Camera
{
    
    // Make all the protected methods that we want to test publically available
    
    public:
        MyCamera(ConfigurationParameters &configParam, HDF5File &hdf5file, Telescope &telescope, Sky &sky)
            : Camera(configParam, hdf5file, telescope, sky) {};

        pair<double, double> test_skyToFocalPlaneCoordinates(double raStar, double decStar)
            {return skyToFocalPlaneCoordinates(raStar, decStar);};
        pair<double, double> test_focalPlaneToSkyCoordinates(double xFPprime, double yFPprime)
            {return focalPlaneToSkyCoordinates(xFPprime, yFPprime);};

        pair<double, double> test_polarToCartesianFocalPlaneCoordinates(double distance, double angle)
            {return polarToCartesianFocalPlaneCoordinates(distance, angle);};
        pair<double, double> test_cartesianToPolarFocalPlaneCoordinates(double xFPdist, double yFPdist)
            {return distortedToUndistortedFocalPlaneCoordinates(xFPdist, yFPdist);};

        pair<double, double> test_undistortedToDistortedFocalPlaneCoordinates(double xFPmm, double yFPmm)
            {return undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm);};
        pair<double, double> test_distortedToUndistortedFocalPlaneCoordinates(double xFPdist, double yFPdist)
            {return distortedToUndistortedFocalPlaneCoordinates(xFPdist, yFPdist);};

        double test_getGnomonicRadialDistanceFromOpticalAxis(double xFPprime, double yFPprime)
            {return getGnomonicRadialDistanceFromOpticalAxis(xFPprime, yFPprime);};

        void test_setDistortionPolynomial(Polynomial1D &polynomial, Polynomial1D &inversePolynomial)
            {setDistortionPolynomial(polynomial, inversePolynomial);};
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

    JitterGenerator *jitterGenerator = new JitterFromRedNoise(configParams);
    DriftGenerator *driftGenerator = new ThermoElasticDriftFromRedNoise(configParams);
    
    Platform platform = Platform(configParams, hdf5File, *jitterGenerator);
    Telescope telescope = Telescope(configParams, hdf5File, platform, *driftGenerator);
    Sky sky = Sky(configParams);

    MyCamera camera = MyCamera(configParams, hdf5File, telescope, sky);

    double raStar, decStar;      // [rad]
    double xFPprime, yFPprime;   // [mm]
    double xFP, yFP;             // [mm]
    double xDeg, yDeg;           // [deg]
    double xFPrad, yFPrad;       // [rad]
    double radius;               // [rad]

    for (auto &data: gnomonic)
    {
        // Other tests that need to be performed is the conversion from xFP, yFP to xDeg, yDeg (in the data table above)

        xDeg = data["xDeg"];
        yDeg = data["yDeg"];
        xFP  = data["xFP"];
        yFP  = data["yFP"];
        xFPrad = deg2rad(xDeg);
        yFPrad = deg2rad(yDeg);

        // We do not use angular focal plane coordinates anymore, but even then, conversion from xDeg, yDeg to xFP, yFP
        // is still an unknown....

        //tie(xFPprime, yFPprime) = camera.test_angularToPlanarFocalPlaneCoordinates(xFPrad, yFPrad);

        //EXPECT_NEAR(xFP, xFPprime, 0.00001);
        //EXPECT_NEAR(yFP, yFPprime, 0.00001);

        //Log.debug("CameraTest.GnomonicRadialDistance: xFPprime, yFPprime [mm]= " + dtos(xFPprime) + ", " + dtos(yFPprime));

        radius = camera.test_getGnomonicRadialDistanceFromOpticalAxis(xFP, yFP); // [mm] -> [radians]

        EXPECT_NEAR(data["radius"], rad2deg(radius), 0.0001);

        //Log.debug("CameraTest.GnomonicRadialDistance: xDeg, yDeg [deg] = " + dtos(xDeg) + ", " + dtos(yDeg));
        //Log.debug("CameraTest.GnomonicRadialDistance: xFPrad, yFPrad [rad] = " + dtos(xFPrad) + ", " + dtos(yFPrad));
        //Log.debug("CameraTest.GnomonicRadialDistance: radius [rad] = " + dtos(radius));
        //Log.debug("CameraTest.GnomonicRadialDistance: radius [deg] = " + dtos(rad2deg(radius)));

        tie(raStar, decStar) = camera.test_focalPlaneToSkyCoordinates(xFP, yFP);

        //Log.debug("CameraTest.GnomonicRadialDistance: raStar, decStar [rad] = " + dtos(raStar) + ", " + dtos(decStar));
        //Log.debug("CameraTest.GnomonicRadialDistance: raStar, decStar [deg] = " + dtos(rad2deg(raStar)) + ", " + dtos(rad2deg(decStar)));

    }

}


// Perform some simple tests providing a simple polynomial as distortion function
// to the Camera class.

TEST_F(CameraTest, distortedCoordinates)
{
    LOG_STARTING_OF_TEST

    using StringUtilities::dtos;

    JitterGenerator *jitterGenerator = new JitterFromRedNoise(configParams);
    DriftGenerator *driftGenerator = new ThermoElasticDriftFromRedNoise(configParams);
    Platform platform = Platform(configParams, hdf5File, *jitterGenerator);
    Telescope telescope = Telescope(configParams, hdf5File, platform, *driftGenerator);
    Sky sky = Sky(configParams);


    
    
    // Simple linear function: y = x
    
    int degree = 1;
    vector<double> coeff = {0.0, 1.0};
    Polynomial1D polynomial = Polynomial1D(degree, coeff);

    MyCamera camera = MyCamera(configParams, hdf5File, telescope, sky);

    camera.test_setDistortionPolynomial(polynomial, polynomial); // Do not care about the inverse polynomial for this test

    double xFPdist, yFPdist;   // [mm]

    tie(xFPdist, yFPdist) = camera.test_undistortedToDistortedFocalPlaneCoordinates(10.0, 0.0);
    EXPECT_NEAR(10.0000, xFPdist, 0.00001);
    EXPECT_NEAR( 0.0000, yFPdist, 0.00001);

    tie(xFPdist, yFPdist) = camera.test_undistortedToDistortedFocalPlaneCoordinates(0.0, 10.0);
    EXPECT_NEAR( 0.0000, xFPdist, 0.00001);
    EXPECT_NEAR(10.0000, yFPdist, 0.00001);

    tie(xFPdist, yFPdist) = camera.test_undistortedToDistortedFocalPlaneCoordinates(5.0, 5.0);
    EXPECT_NEAR( 5.0000, xFPdist, 0.00001);
    EXPECT_NEAR( 5.0000, yFPdist, 0.00001);


    
    
    // Polynomial function: y = 2.0 + 0.5 * x + 1.5 * x^2
    
    degree = 2;
    coeff = {2.0, 0.5, 1.5};
    polynomial = Polynomial1D(degree, coeff);

    camera.test_setDistortionPolynomial(polynomial, polynomial); // Do not care about the inverse polynomial for this test

    tie(xFPdist, yFPdist) = camera.test_undistortedToDistortedFocalPlaneCoordinates(10.0, 0.0);
    EXPECT_NEAR(157.0000, xFPdist, 0.00001);
    EXPECT_NEAR(  0.0000, yFPdist, 0.00001);

    tie(xFPdist, yFPdist) = camera.test_undistortedToDistortedFocalPlaneCoordinates(0.0, 10.0);
    EXPECT_NEAR(  0.0000, xFPdist, 0.00001);
    EXPECT_NEAR(157.0000, yFPdist, 0.00001);

    tie(xFPdist, yFPdist) = camera.test_undistortedToDistortedFocalPlaneCoordinates(5.0, 5.0);
    EXPECT_NEAR( 56.947222, xFPdist, 0.00001);
    EXPECT_NEAR( 56.947222, yFPdist, 0.00001);

}


TEST_F(CameraTest, reproduceDistortionMap)
{
    LOG_STARTING_OF_TEST

    using StringUtilities::dtos;

    JitterGenerator *jitterGenerator = new JitterFromRedNoise(configParams);
    DriftGenerator *driftGenerator = new ThermoElasticDriftFromRedNoise(configParams);
    Platform platform = Platform(configParams, hdf5File, *jitterGenerator);
    Telescope telescope = Telescope(configParams, hdf5File, platform, *driftGenerator);
    Sky sky = Sky(configParams);

    // Just a few values from the current distortion map for which the Polynomial1D was fitted.
    
    vector<map<string, double>> distortion;
    distortion.push_back(map<string, double> {{"xFPmm",  0.000000}, {"yFPmm",  0.000000}, {"xFPdist",   0.000000}, {"yFPdist",  0.000000}});

    distortion.push_back(map<string, double> {{"xFPmm",  2.156636}, {"yFPmm", 12.951322}, {"xFPdist",   2.158552}, {"yFPdist", 12.962833}});  // 0.500000	3.000000	3.041231
    distortion.push_back(map<string, double> {{"xFPmm", 12.951322}, {"yFPmm",  2.156636}, {"xFPdist",  12.962833}, {"yFPdist",  2.158552}});  // 3.000000	0.500000	3.041231

    distortion.push_back(map<string, double> {{"xFPmm", 22.490223}, {"yFPmm", 52.077606}, {"xFPdist",  22.869180}, {"yFPdist", 52.955106}});  // 5.200000   11.900000   12.927980
    distortion.push_back(map<string, double> {{"xFPmm", 52.077606}, {"yFPmm", 22.490223}, {"xFPdist",  52.955106}, {"yFPdist", 22.869180}});  // 11.900000  5.200000    12.927980

    distortion.push_back(map<string, double> {{"xFPmm", 10.789752}, {"yFPmm",  0.000000}, {"xFPdist",  10.796226}, {"yFPdist",  0.000000}});
    distortion.push_back(map<string, double> {{"xFPmm", 24.666449}, {"yFPmm",  0.000000}, {"xFPdist",  24.743989}, {"yFPdist",  0.000000}});
    distortion.push_back(map<string, double> {{"xFPmm", 68.998805}, {"yFPmm",  0.000000}, {"xFPdist",  70.734282}, {"yFPdist",  0.000000}});
    distortion.push_back(map<string, double> {{"xFPmm", 80.296089}, {"yFPmm",  0.000000}, {"xFPdist",  83.062336}, {"yFPdist",  0.000000}});


    
    
    
    // These values are for a fit of Polynomial1D with degree=3 to the distortion table

    int degree = 3;
    vector<double> coeff = {-0.0036696919678, 1.0008542317, -4.12553764817e-05, 5.7201219949e-06};
    Polynomial1D polynomial = Polynomial1D(degree, coeff);

    vector<double> inverseCoeff = {-0.00458067036444, 1.00110311283, -5.61136295937e-05, -4.311925329e-06};
    Polynomial1D inversePolynomial = Polynomial1D(degree, inverseCoeff);

    MyCamera camera = MyCamera(configParams, hdf5File, telescope, sky);
    camera.test_setDistortionPolynomial(polynomial, inversePolynomial);

    double xFPdist, yFPdist;   // [mm]
    double xFPmm, yFPmm;       // [mm]

    for (auto &data: distortion)
    {
        xFPmm = data["xFPmm"];
        yFPmm = data["yFPmm"];

        tie(xFPdist, yFPdist) = camera.test_undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm);
        EXPECT_NEAR( data["xFPdist"], xFPdist, 0.006);
        EXPECT_NEAR( data["yFPdist"], yFPdist, 0.006);

        xFPdist = data["xFPdist"];
        yFPdist = data["yFPdist"];

        tie(xFPmm, yFPmm) = camera.test_distortedToUndistortedFocalPlaneCoordinates(xFPdist, yFPdist);
        EXPECT_NEAR( data["xFPmm"], xFPmm, 0.006);
        EXPECT_NEAR( data["yFPmm"], yFPmm, 0.006);

    }

    
    // The following loop is to test if the inverse function is indeed the inverse polynomial.
    // The biggest difference that I got was 0.0103 mm which is only slightly less than the size of a pixel.
    // TODO: Should we look into improving this distortion with Polynomials?

    for (auto &data: distortion)
    {
        double xFPmm = data["xFPmm"];
        double yFPmm = data["yFPmm"];

        double xFPdist = polynomial(xFPmm);
        double yFPdist = polynomial(yFPmm);

        EXPECT_NEAR(xFPmm, inversePolynomial(xFPdist), 0.02);
        EXPECT_NEAR(yFPmm, inversePolynomial(yFPdist), 0.02);
    }


    
    

    
    // These values are for a fit of Polynomial1D with degree=4 to the distortion table.
    // Obviously the error of the fit is smaller, but we need to check what the effect is on the noise level.
    // Remember that 1 pixel = 0.018 mm

    degree = 4;
    coeff = {0.000814670391532, 0.99970324817, 2.39592182367e-05, 4.44973838376e-06, 7.93413878401e-09};
    polynomial = Polynomial1D(degree, coeff);

    camera.test_setDistortionPolynomial(polynomial, polynomial); // Do not care about the inverse polynomial for this test

    for (auto &data: distortion)
    {
        double xFPmm = data["xFPmm"];
        double yFPmm = data["yFPmm"];

        tie(xFPdist, yFPdist) = camera.test_undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm);
        EXPECT_NEAR( data["xFPdist"], xFPdist, 0.002);
        EXPECT_NEAR( data["yFPdist"], yFPdist, 0.002);
    }

}




