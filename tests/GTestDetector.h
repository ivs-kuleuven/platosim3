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
	void test_setSubPixelMap(const arma::fmat &subPixelMap)
	{
		if ((subPixelMap.n_rows != this->subPixelMap.n_rows)
				|| (subPixelMap.n_cols != this->subPixelMap.n_cols))
		{
			Log.error(
					"MyDetector: test_setSubPixelMap with incompatible array shape");
			exit(1);
		}

		this->subPixelMap = subPixelMap;
	}
	;

	void test_setBiasRegisterMap(const arma::fmat &biasMap)
	{
		if ((biasMap.n_rows != this->biasMap.n_rows)
				|| (biasMap.n_cols != this->biasMap.n_cols))
		{
			Log.error(
					"MyDetector: test_setBiasRegisterMap with incompatible array shape");
			exit(1);
		}

		this->biasMap = biasMap;
	}
	;
	void test_setSmearingMap(const arma::fmat &smearingMap)
	{
		if ((smearingMap.n_rows != this->smearingMap.n_rows)
				|| (smearingMap.n_cols != this->smearingMap.n_cols))
		{
			Log.error(
					"MyDetector: test_setSmearingMap with incompatible array shape");
			exit(1);
		}

		this->smearingMap = smearingMap;
	};

	arma::Mat<float> test_getSubfield()
	{
		return getSubfield();
	}
	;

	arma::fmat test_getSubPixelMap()
	{
		return subPixelMap;
	};

	arma::fmat test_getBiasRegisterMap()
	{
		return biasMap;
	};

	arma::fmat test_getSmearingMap()
	{
		return smearingMap;
	};
	arma::fmat test_getFlatfieldMap()
	{
		return flatfieldMap;
	};

	void test_reset(){
		reset();
	};

	void test_addElectronicOffset()
	{
		addElectronicOffset();
	};

	void test_applyGain()
	{
		applyGain();
	};

	void test_applyQuantumEfficiency()
	{
		applyQuantumEfficiency();
	}

	void test_generateFlatFieldMap()
	{
		generateFlatfieldMap();
	}

	void test_applyFlatfield()
	{
		applyFlatfield();
	}

	bool test_isInSubPixelMap(double row, double column)
	{
		return isInSubPixelMap(row, column);
	}
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










/**
 * Dimensions.
 *
 * The pixel map, bias register map, and smearing map must be generated at pixel level, whilst the
 * sub-pixel map and the flatfield map must be generated at sub-pixel level.
 *
 * Currently, no edge pixels are added to the sub-pixel maps.  At a later stage, the addition of edge pixels
 * will be implemented and the test harness must be updated accordingly.
 */
TEST_F(DetectorTest, dimensions)
{
	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");
	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	// Sub-pixel map
	// TODO Should initially also include edge pixels (not implemented currently)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	// Pixel map

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

	// Bias register map

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	// Smearing map

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

	// Flatfield map

	ASSERT_DOUBLE_EQ(numRowsSubField * numSubPixels, detector.test_getFlatfieldMap().n_rows);
	ASSERT_DOUBLE_EQ(numColumnsSubField * numSubPixels, detector.test_getFlatfieldMap().n_cols);
}










TEST_F(DetectorTest, generateFlatfield)
{

}










/**
 * Reset.
 *
 * The dimensions of the sub-pixel map, pixel map, bias register map, and smearing map must remain unchanged
 * but the values in these maps must be set to zero.
 */
TEST_F(DetectorTest, reset)
{
	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");
	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	detector.test_setSubPixelMap(arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels));
	detector.test_setSubfield(arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField));
	detector.test_setBiasRegisterMap(arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField));
	detector.test_setSmearingMap(arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField));

	detector.test_reset();

	// Sub-pixel map

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);
	ASSERT_EQ(0, arma::accu(arma::abs(detector.test_getSubPixelMap())));

	// Pixel map

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);
	ASSERT_EQ(0, arma::accu(arma::abs(detector.test_getSubfield())));

	// Bias register map

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);
	ASSERT_EQ(0, arma::accu(arma::abs(detector.test_getBiasRegisterMap())));

	// Smearing map

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);
	ASSERT_EQ(0, arma::accu(arma::abs(detector.test_getSmearingMap())));

	// Flatfield map not reset!

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getFlatfieldMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getFlatfieldMap().n_cols);
	ASSERT_NE(0, arma::accu(arma::abs(detector.test_getFlatfieldMap())));
}









/**
 * Flatfielding.
 *
 * The sub-pixel map must be divided by the flatfield map.
 */
TEST_F(DetectorTest, applyFlatfield)
{
	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");
	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const double quantumEfficiency = configParams.getDouble("CCD/QuantumEfficiency");

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);

	detector.test_generateFlatFieldMap();
	arma::fmat flatfieldMap = detector.test_getFlatfieldMap();
	detector.test_applyFlatfield();

	// Sub-pixel map
	// TODO Should also include edge pixels (not implemented currently)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSubPixelMap()) == arma::vectorise(subPixelMap / flatfieldMap)));

	// Pixel map

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSubfield()) == arma::vectorise(subField)));

	// Bias register map

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getBiasRegisterMap()) == arma::vectorise(biasMap)));

	// Smearing map

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSmearingMap()) == arma::vectorise(smearingMap)));
}









TEST_F(DetectorTest, rebin)
{

}








TEST_F(DetectorTest, addFlux)
{
	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);


}









TEST_F(DetectorTest, isInSubField)
{

}









/**
 * Check whether given position (row, column) is located in the sub-pixel map.
 */
TEST_F(DetectorTest, isInSubPixelMap)
{
	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	ASSERT_FALSE(detector.test_isInSubPixelMap(-1,0));
	ASSERT_FALSE(detector.test_isInSubPixelMap(-1,-20));
	ASSERT_FALSE(detector.test_isInSubPixelMap(-1,-50));
	ASSERT_FALSE(detector.test_isInSubPixelMap(-1,0));

	ASSERT_TRUE(detector.test_isInSubPixelMap(numRowsSubField * numSubPixels / 2.0, numColumnsSubField * numSubPixels / 2.0));

	ASSERT_FALSE(detector.test_isInSubPixelMap(numRowsSubField * numSubPixels / 2.0, numColumnsSubField * numSubPixels));
	ASSERT_FALSE(detector.test_isInSubPixelMap(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels / 2.0));

	ASSERT_TRUE(detector.test_isInSubPixelMap(numRowsSubField * numSubPixels -1.0, numColumnsSubField * numSubPixels - 1.0));

	ASSERT_FALSE(detector.test_isInSubPixelMap(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels - 1.0));
	ASSERT_FALSE(detector.test_isInSubPixelMap(numRowsSubField * numSubPixels -1.0, numColumnsSubField * numSubPixels));
}









/**
 * Quantum efficiency.
 *
 * Multiplies the values in the sub-pixel map with the quantum efficiency.
 */
TEST_F(DetectorTest, applyQuantumEfficiency)
{
	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");
	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const double quantumEfficiency = configParams.getDouble("CCD/QuantumEfficiency");

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);

	detector.test_applyQuantumEfficiency();

	// Sub-pixel map

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSubPixelMap()) == arma::vectorise(subPixelMap)));

	// Pixel map

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSubfield()) == arma::vectorise(subField * quantumEfficiency)));

	// Bias register map

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getBiasRegisterMap()) == arma::vectorise(biasMap)));

	// Smearing map

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSmearingMap()) == arma::vectorise(smearingMap)));
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









/**
 * Gain.
 *
 * The values in the pixel map, bias register map, and smearing map must be divided by the
 * gain, on order to convert from [e- / pixel] to [ADU / pixel].
 */
TEST_F(DetectorTest, applyGain)
{
	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int gain = configParams.getInteger("CCD/Gain");

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);

	detector.test_applyGain();

	// Sub-pixel map

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSubPixelMap()) == arma::vectorise(subPixelMap)));

	// Pixel map

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSubfield()) == arma::vectorise(subField / gain)));

	// Bias register map

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getBiasRegisterMap()) == arma::vectorise(biasMap / gain)));

	// Smearing map

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSmearingMap()) == arma::vectorise(smearingMap / gain)));
}









/**
 * Electronic offset.
 *
 * The electronic offset must be added to the pixel map, bias register map, and smearing map.
 */
TEST_F(DetectorTest, addElectronicOffset)
{
	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int electronicOffset = configParams.getInteger("CCD/ElectronicOffset");

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);

	detector.test_addElectronicOffset();

	// Sub-pixel map

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSubPixelMap()) == arma::vectorise(subPixelMap)));

	// Pixel map

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSubfield()) == arma::vectorise(subField + electronicOffset)));

	// Bias register map

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getBiasRegisterMap()) == arma::vectorise(biasMap + electronicOffset)));

	// Smearing map

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);
	EXPECT_TRUE(arma::all(arma::vectorise(detector.test_getSmearingMap()) == arma::vectorise(smearingMap + electronicOffset)));
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

