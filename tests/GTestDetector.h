#include <cstdio>
#include <cmath>
#include <map>

#include "armadillo"

#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "Simulation.h"


using namespace std;



/**
 * \class DetectorTest
 * 
 * \brief Test Fixture for the Detector class. 
 * 
 * \details
 * 
 * Setup configuration parameters for the Detector class and handle the HDF5 open-close.
 * 
 * Input parameters for this test are located in the file 'input_DetectorTest.yaml' in the
 * testData directory of the distribution.
 * 
 * The test creates an HDF5 file 'detectorTest.hdf5' in the current working directory. This file is removed after the test finishes.
 */

class DetectorTest: public testing::Test
{
protected:

	virtual void SetUp()
	{
		configParams = ConfigurationParameters(
				"../testData/input_DetectorTest.yaml");
		hdf5File.open(hdf5Filename);
	}

	virtual void TearDown()
	{
		hdf5File.close();
		FileUtilities::remove(hdf5Filename);
	}

	void reset()
	{
		hdf5File.close();
		FileUtilities::remove(hdf5Filename);
		hdf5File.open(hdf5Filename);
	}

	string hdf5Filename = "detectorTest.hdf5";
	ConfigurationParameters configParams;
	HDF5File hdf5File;
};












/**
 * 
 * \brief This subclass of Detector serves the sole purpose of testing protected methods of Detector.
 * 
 */

class MyDetector: public Detector
{
public:
	MyDetector(ConfigurationParameters &configParam, HDF5File &hdf5File,
			Camera &camera) :
			Detector(configParam, hdf5File, camera)
	{
	}
	;

	pair<double, double> test_pixelToPlanarFocalPlaneCoordinates(double row,
			double column)
	{
		return pixelToPlanarFocalPlaneCoordinates(row, column);
	}
	;
	pair<double, double> test_planarFocalPlaneToPixelCoordinates(
			double xFPprime, double yFPprime)
	{
		return planarFocalPlaneToPixelCoordinates(xFPprime, yFPprime);
	}
	;

	void test_setSubfield(const arma::Mat<float> &subfield)
	{
		setSubfield(subfield);
	}
	;
	arma::Mat<float> test_getSubfield()
	{
		return getSubfield();
	}
	;

};












TEST_F(DetectorTest, checkConversionsBetweenPixelsAndFocalPlane)
{

    LOG_STARTING_OF_TEST

    double row, column;
    double xFPprime, yFPprime;

    vector<map<string, double>> pixel2fp;

    // Initialise all objects necessary to set up a Detector object

    JitterFromRedNoise jitterGenerator(configParams);
    Sky sky(configParams);

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
        configParams.setParameter("CCD/OriginOffsetX", to_string(data["zeroPointX"]));
        configParams.setParameter("CCD/OriginOffsetY", to_string(data["zeroPointY"]));
        configParams.setParameter("CCD/Orientation", to_string(data["ccdAngle"]));
    
        Platform platform(configParams, hdf5File, jitterGenerator);
        Telescope telescope(configParams, hdf5File, platform);
        Camera camera(configParams, hdf5File, telescope, sky);
        MyDetector detector(configParams, hdf5File, camera);
    
        row = data["xCCD"];
        column = data["yCCD"];
        tie(xFPprime, yFPprime) = detector.test_pixelToPlanarFocalPlaneCoordinates(row, column);
    
        EXPECT_NEAR(data["xFP"], xFPprime, 0.00001);
        EXPECT_NEAR(data["yFP"], yFPprime, 0.00001);    
    }
}











TEST_F(DetectorTest, setAndGetSubfield)
{
	LOG_STARTING_OF_TEST

	// Initialise all objects necessary to set up a Detector object

JitterFromRedNoise	jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Find out what size of subfield we specified in the input yaml file

	const int Nrows = configParams.getInteger("SubField/NumRows");
	const int Ncols = configParams.getInteger("SubField/NumColumns");

	// Make our own subfield (a diagonal unity matrix) and feed it to detector

	auto diagonalMatrix = arma::eye<arma::Mat<float>>(Nrows, Ncols);
	detector.test_setSubfield(diagonalMatrix);

	// Get back the subfield from the detector

	arma::Mat<float> mySubfield = detector.test_getSubfield();

	// Compare whether the output is the same matrix as we put in

	EXPECT_EQ(mySubfield.n_rows, Nrows);
	EXPECT_EQ(mySubfield.n_cols, Ncols);
	EXPECT_TRUE(arma::all(arma::vectorise(mySubfield) == arma::vectorise(diagonalMatrix)));

}










TEST_F(DetectorTest, generateFlatfield)
{

}










TEST_F(DetectorTest, reset)
{

}









TEST_F(DetectorTest, applyFlatfield)
{

}









TEST_F(DetectorTest, rebin)
{

}








TEST_F(DetectorTest, addFlux)
{

}









TEST_F(DetectorTest, isInSubField)
{

}









TEST_F(DetectorTest, isInSubPixelMap)
{

}









TEST_F(DetectorTest, applyQuantumEfficiency)
{

}









TEST_F(DetectorTest, addPhotonNoise)
{

}









TEST_F(DetectorTest, applyFullWellSaturation)
{

}









TEST_F(DetectorTest, applyCte)
{

}









TEST_F(DetectorTest, applyOpenShutterSmearing)
{

}









TEST_F(DetectorTest, addReadoutNoise)
{

}









TEST_F(DetectorTest, applyGain)
{

}









TEST_F(DetectorTest, addElectronicOffset)
{

}









TEST_F(DetectorTest, applyDigitalSaturation)
{

}









TEST_F(DetectorTest, convolveWithPsf)
{

}









TEST_F(DetectorTest, getPlanarFocalPlaneCoordinatesOfSubfieldCorners)
{

}









TEST_F(DetectorTest, getSolidAngleOfOnePixel)
{

}

