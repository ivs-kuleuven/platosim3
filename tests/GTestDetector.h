#include <cstdio>
#include <cmath>
#include <map>
#include <random>
#include <algorithm>

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

	arma::fmat applyCteOldImplementation(arma::fmat pixelMap, int subFieldZeroPointRow, int numRowsPixelMap, int numColumnsPixelMap, double meanCte)
	{
		// Create a map in which we will shift the rows of the pixel map one-by-one
		// towards the readout register.  Bear in mind that the bottom row of the
		// sub-field is not necessarily right next to the readout register (the
		// distance between the two is subFieldZeroPointRow).

		arma::Mat<float> shiftMap;
		shiftMap.zeros(subFieldZeroPointRow + numRowsPixelMap, numColumnsPixelMap);
		shiftMap.submat(arma::span(subFieldZeroPointRow, subFieldZeroPointRow + numRowsPixelMap - 1), arma::span::all) = pixelMap;

		arma::fmat cteMap (numRowsPixelMap + subFieldZeroPointRow, numColumnsPixelMap);
		cteMap.fill(meanCte);

		// The readout register

		arma::Row<float> readoutStrip;
		readoutStrip.zeros(numColumnsPixelMap);

		// Array filled with ones (needed for the CTI)

		arma::Row<float> ones;
		ones.ones(numColumnsPixelMap);

		// Shift all the rows down (i.e. towards the readout register) one-by-one
		// Keep on doing this until all rows have been read out.

		for (int shiftIndex = 0; shiftIndex < numRowsPixelMap + subFieldZeroPointRow; shiftIndex++)
		{

			// Shift the bottom row to the readout strip

			readoutStrip = cteMap(0, arma::span::all) % shiftMap(0, arma::span::all);

			if (shiftIndex >= subFieldZeroPointRow)
			{
				pixelMap(shiftIndex - subFieldZeroPointRow, arma::span::all) = readoutStrip(0, arma::span::all);
			}

			// Shift all other rows one row down (i.e. closer to the readout register)

			for (int row = 0; row < numRowsPixelMap + subFieldZeroPointRow - 1; row++)
			{
				shiftMap(row, arma::span::all) = (ones
						- cteMap(row, arma::span::all))
						% shiftMap(row, arma::span::all)// Left behind when shifting row down (CTI = 1 - CTE)
						+ cteMap(row + 1, arma::span::all)
								% shiftMap(row + 1, arma::span::all);// Transferred (CTE)
				}
			}

		return pixelMap;
	}
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



	pair<double, double> test_pixelToPlanarFocalPlaneCoordinates(double row,
			double column)
	{
		return pixelToPlanarFocalPlaneCoordinates(row, column);
	}



	pair<double, double> test_planarFocalPlaneToPixelCoordinates(
			double xFPprime, double yFPprime)
	{
		return planarFocalPlaneToPixelCoordinates(xFPprime, yFPprime);
	}



	void test_setSubfield(const arma::Mat<float> &subfield)
	{
		setSubfield(subfield);
	}



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
	}



	arma::Mat<float> test_getSubfield()
	{
		return getSubfield();
	}



	arma::fmat test_getSubPixelMap()
	{
		return subPixelMap;
	}



	arma::fmat test_getBiasRegisterMap()
	{
		return biasMap;
	}



	arma::fmat test_getSmearingMap()
	{
		return smearingMap;
	}



	arma::fmat test_getFlatfieldMap()
	{
		return flatfieldMap;
	}



	void test_reset(){
		reset();
	}



	void test_addElectronicOffset()
	{
		addElectronicOffset();
	}



	void test_applyGain()
	{
		applyGain();
	}



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



	void test_rebin()
	{
		rebin();
	}



	void test_applyDigitalSaturation()
	{
		applyDigitalSaturation();
	}



	void test_applyCte()
	{
		applyCte();
	}



	void test_addPhotonNoise()
	{
		addPhotonNoise();
	}

	void test_addReadoutNoise()
	{
		addReadoutNoise();
	}



	void test_applyOpenShutterSmearing(double exposureTime)
	{
		applyOpenShutterSmearing(exposureTime);
	}



	void test_setZeroPointRow(int row)
	{
		subFieldZeroPointRow = row;
	}



	normal_distribution<double> test_getReadoutNoiseDistribution()
	{
		return readoutNoiseDistribution;
	}



	void test_setDimensions(unsigned int numRowsSubField, unsigned int numColumnsSubField, unsigned numSubPixels)
	{
		numRowsPixelMap = numRowsSubField;
		numColumnsPixelMap = numColumnsSubField;

		numRowsSubPixelMap = numRowsSubField * numSubPixels;
		numColumnsSubPixelMap = numColumnsSubField * numSubPixels;
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
        // FIXME: This is an inconvenience that the selection/orientation of the CCD can not be done with access methods.
        
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
	LOG_STARTING_OF_TEST

	// Constructiomn

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");



	// Sub-pixel map: check dimensions

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	// Pixel map: check dimensions

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

	// Bias register map: check dimensions

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	// Smearing map: check dimensions

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

	// Flatfield map: check dimensions

	ASSERT_DOUBLE_EQ(numRowsSubField * numSubPixels, detector.test_getFlatfieldMap().n_rows);
	ASSERT_DOUBLE_EQ(numColumnsSubField * numSubPixels, detector.test_getFlatfieldMap().n_cols);
}










TEST_F(DetectorTest, generateFlatfield)
{
	LOG_STARTING_OF_TEST

}










/**
 * Reset.
 *
 * The dimensions of the sub-pixel map, pixel map, bias register map, and smearing map must remain unchanged
 * but the values in these maps must be set to zero.
 */
TEST_F(DetectorTest, reset)
{
	// Construction

	LOG_STARTING_OF_TEST

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	// Initialise sub-pixel map, pixel map, bias register map, and smearing map

	detector.test_setSubPixelMap(arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels));
	detector.test_setSubfield(arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField));
	detector.test_setBiasRegisterMap(arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField));
	detector.test_setSmearingMap(arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField));



	// Reset

	detector.test_reset();



	// Sub-pixel map: check dimensions and content (should be all zeroes)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	ASSERT_EQ(0, arma::accu(arma::abs(detector.test_getSubPixelMap())));

	// Pixel map: check dimensions and content (should be all zeroes)

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

	ASSERT_EQ(0, arma::accu(arma::abs(detector.test_getSubfield())));

	// Bias register map: check dimensions and content (should be all zeroes)

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	ASSERT_EQ(0, arma::accu(arma::abs(detector.test_getBiasRegisterMap())));

	// Smearing map: check dimensions and content (should be all zeroes)

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
	LOG_STARTING_OF_TEST

	// Construction

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	// Initialise sub-pixel map, pixel map, bias register map, and smearing map

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);
	EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap) == arma::vectorise(detector.test_getSubPixelMap())));
	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);



	// Flatfield

	detector.test_applyFlatfield();



	// Sub-pixel map: check dimensions and content (divided by flatfield map)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap % detector.test_getFlatfieldMap()) == arma::vectorise(detector.test_getSubPixelMap())));

	// Pixel map: check dimensions and content (unaltered)

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subField) == arma::vectorise(detector.test_getSubfield())));

	// Bias register map: check dimensions and content (unaltered)

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(biasMap) == arma::vectorise(detector.test_getBiasRegisterMap())));

	// Smearing map: check dimensions and content (unaltered)

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(smearingMap) == arma::vectorise(detector.test_getSmearingMap())));
}









/**
 * Rebinning.
 */
TEST_F(DetectorTest, rebin)
{
	LOG_STARTING_OF_TEST

	// Construction

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	// Initialise sub-pixel map, pixel map, bias register map, and smearing map

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	subPixelMap.ones(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);



	// Rebin

	detector.test_rebin();



	// Pixel map: check dimensions  and content (all sub-pixels must be summed per pixel)

	ASSERT_EQ(detector.test_getSubPixelMap().n_rows / numSubPixels, detector.test_getSubfield().n_rows);
	ASSERT_EQ(detector.test_getSubPixelMap().n_cols / numSubPixels, detector.test_getSubfield().n_cols);

	double expectedValue;

	for (unsigned int row = 0; row < numRowsSubField; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			expectedValue = 0.0;

			for(unsigned int subRow = row * numSubPixels; subRow < (row + 1) * numSubPixels; subRow++)
			{
				for(unsigned int subColumn = column * numSubPixels; subColumn < (column + 1) * numSubPixels; subColumn++)
				{
					expectedValue += subPixelMap(row, column);
				}
			}

			ASSERT_EQ(expectedValue, detector.test_getSubfield()(row, column));
		}
	}

	// Sub-pixel map: check dimensions and content (unaltered)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap) == arma::vectorise(detector.test_getSubPixelMap())));

	// Bias register map: check dimensions and content (unaltered)

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(biasMap) == arma::vectorise(detector.test_getBiasRegisterMap())));

	// Smearing map: check dimensions and content (unaltered)

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(smearingMap) == arma::vectorise(detector.test_getSmearingMap())));
}








/**
 * Adding same flux value to all sub-pixels.
 */
TEST_F(DetectorTest, addBackgroudFlux)
{
	LOG_STARTING_OF_TEST

	// Construction

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	// Initialise sub-pixel map, pixel map, bias register map, and smearing map

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);



	// Add flux

	double background = 121.47;	// Per pixel

	detector.addFlux(background);



	// Sub-pixel map: check dimensions and content (background must be added)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap + (background / numSubPixels / numSubPixels)) == arma::vectorise(detector.test_getSubPixelMap())));

	// Pixel map: check dimensions and content (unaltered)

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subField) == arma::vectorise(detector.test_getSubfield())));

	// Bias register map: check dimensions and content (unaltered)

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(biasMap) == arma::vectorise(detector.test_getBiasRegisterMap())));

	// Smearing map: check dimensions and content (unaltered)

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(smearingMap) == arma::vectorise(detector.test_getSmearingMap())));
}









TEST_F(DetectorTest, addFlux)
{
	LOG_STARTING_OF_TEST

}









TEST_F(DetectorTest, isInSubField)
{
	LOG_STARTING_OF_TEST

}









/**
 * Check whether given position (row, column) is located in the sub-pixel map.
 */
TEST_F(DetectorTest, isInSubPixelMap)
{
	LOG_STARTING_OF_TEST

	// Construction

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

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
	LOG_STARTING_OF_TEST

	// Construction

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	const double quantumEfficiency = configParams.getDouble("CCD/QuantumEfficiency");

	// Initialise sub-pixel map, pixel map, bias register map, and smearing map

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);



	// Quantum efficiency

	detector.test_applyQuantumEfficiency();



	// Sub-pixel map: check dimensions and content (unaltered)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap) == arma::vectorise(detector.test_getSubPixelMap())));

	// Pixel map: check dimensions and content (multiplied by quantum efficiency)

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subField * quantumEfficiency) == arma::vectorise(detector.test_getSubfield())));

	// Bias register map: check dimensions and content (unaltered)

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(biasMap) == arma::vectorise(detector.test_getBiasRegisterMap())));

	// Smearing map: check dimensions and content (unaltered)

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(smearingMap) == arma::vectorise(detector.test_getSmearingMap())));
}









/**
 * Photon noise.
 *
 * Photon noise must be added to the pixel map and the smearing map.
 *
 * As each pixel is treated independently, we repeat the process of adding photon noise (each time to the
 * original pixel map and smearing map) and check afterwards whether this follows the expected Poisson
 * distribution.  We use the normal approximation to the Poisson distribution for testing.
 */
TEST_F(DetectorTest, DISABLED_addPhotonNoise)
{
	LOG_STARTING_OF_TEST

	// Construction

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	const bool includePhotonNoise = configParams.getBoolean("CCD/IncludePhotonNoise");

	// Initialise sub-pixel map, pixel map, bias register map, and smearing map

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::abs(arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField));
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);



	if(includePhotonNoise)
	{
		arma::fmat residualSubField(numRowsSubField, numColumnsSubField);
		arma::fmat meanSubField(numRowsSubField, numColumnsSubField);
		arma::fmat stdDevSubField(numRowsSubField, numColumnsSubField);

		arma::fmat residualSmearingMap(numSmearingOverScanRows, numColumnsSubField);
		arma::fmat meanSmearingMap(numSmearingOverScanRows, numColumnsSubField);
		arma::fmat stdDevSmearingMap(numSmearingOverScanRows, numColumnsSubField);

		int numIterations = 1000;

		for(unsigned int iteration = 0; iteration < numIterations; iteration++)
		{
			// Photon noise

			detector.test_addPhotonNoise();




			residualSubField = detector.test_getSubfield() - subField;
			meanSubField += residualSubField;
			stdDevSubField += (residualSubField % residualSubField);

			residualSmearingMap = detector.test_getSmearingMap() - smearingMap;
			meanSmearingMap += residualSmearingMap;
			stdDevSmearingMap += (residualSmearingMap % residualSmearingMap);

			detector.test_setSubfield(subField);
			detector.test_setSmearingMap(smearingMap);
		}

		meanSubField /= numIterations;
		stdDevSubField /= numIterations;
		stdDevSubField = sqrt(stdDevSubField);

		meanSmearingMap /= numIterations;
		stdDevSmearingMap /= numIterations;
		stdDevSmearingMap = sqrt(stdDevSmearingMap);

		// Pixel map: check dimensions and content (added Poisson distribution, all pixels treated independently)

		EXPECT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
		EXPECT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

		for(unsigned int row = 0; row < numRowsSubField; row++)
		{
			for(unsigned int column = 0; column < numColumnsSubField; column++)
			{
				EXPECT_FLOAT_EQ(sqrt(subField(row, column)), stdDevSubField(row, column));	// Std.dev. = SQRT(pixel value)
			}
		}

		EXPECT_FLOAT_EQ(0.0, meanSubField.min());	// Mean = 0
		EXPECT_FLOAT_EQ(0.0, meanSubField.max());	// Mean = 0

		// Smearing map: check dimensions and content (added Poisson distribution, all pixels treated independently)

		ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
		ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

		for(unsigned int row = 0; row < numSmearingOverScanRows; row++)
		{
			for(unsigned int column = 0; column < numColumnsSubField; column++)
			{
				EXPECT_FLOAT_EQ(sqrt(smearingMap(row, column)), stdDevSmearingMap(row, column));	// Std.dev. = SQRT(pixel value)
			}
		}

		EXPECT_FLOAT_EQ(0.0, meanSmearingMap.min());	// Mean = 0
		EXPECT_FLOAT_EQ(0.0, meanSmearingMap.max());	// Mean = 0


		// Sub-pixel map: check dimensions and content (unaltered)

		ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
		ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

		EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap) == arma::vectorise(detector.test_getSubPixelMap())));


		// Bias register map: check dimensions and content (unaltered)

		ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
		ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

		EXPECT_TRUE(arma::all(arma::vectorise(biasMap) == arma::vectorise(detector.test_getBiasRegisterMap())));

		// Smearing map: check dimensions and content (unaltered)

		EXPECT_TRUE(arma::all(arma::vectorise(smearingMap) == arma::vectorise(detector.test_getSmearingMap())));
	}

	else{

		detector.test_addPhotonNoise();

		// Sub-pixel map: check dimensions and content (unaltered)

		ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
		ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

		EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap) == arma::vectorise(detector.test_getSubPixelMap())));

		// Pixel map: check dimensions and content (unaltered)

		ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
		ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

		EXPECT_TRUE(arma::all(arma::vectorise(subField) == arma::vectorise(detector.test_getSubfield())));

		// Bias register map: check dimensions and content (unaltered)

		ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
		ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

		EXPECT_TRUE(arma::all(arma::vectorise(biasMap) == arma::vectorise(detector.test_getBiasRegisterMap())));

		// Smearing map: check dimensions and content (unaltered)

		ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
		ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

		EXPECT_TRUE(arma::all(arma::vectorise(smearingMap) == arma::vectorise(detector.test_getSmearingMap())));
	}
}









TEST_F(DetectorTest, applyFullWellSaturation)
{
	LOG_STARTING_OF_TEST

}









/**
 * Charge Transfer Efficiency (CTE).
 */
TEST_F(DetectorTest, applyCte)
{
	LOG_STARTING_OF_TEST

	// Construction

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	for(unsigned int detectorZeroPointRow = 0; detectorZeroPointRow <= 15; detectorZeroPointRow += 10)
	{
		detector.test_setZeroPointRow(detectorZeroPointRow);


		// Configuration parameters

		const int numRowsSubField = configParams.getInteger("SubField/NumRows");
		const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

		const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
		const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

		const int numSubPixels = configParams.getInteger("SubField/SubPixels");

		const double meanCte = configParams.getDouble("CCD/CTEMean");

		// Initialise sub-pixel map, pixel map, bias register map, and smearing map

		arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
		detector.test_setSubPixelMap(subPixelMap);

		arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
		//	subField *= 1000.0;
		detector.test_setSubfield(subField);

		arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
		detector.test_setBiasRegisterMap(biasMap);

		arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
		detector.test_setSmearingMap(smearingMap);



		// CTE

		detector.test_applyCte();



		// Sub-pixel map: check dimensions and content (unaltered)

		ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
		ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

		EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap) == arma::vectorise(detector.test_getSubPixelMap())));

		// Pixel map: check dimension and content (compare with brute-force method (i.e. old implementation))

		ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
		ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

		arma::fmat expected = applyCteOldImplementation(subField, detectorZeroPointRow, numRowsSubField, numColumnsSubField, meanCte);

		for(unsigned int row = 0; row < numRowsSubField; row++)
		{
			for(unsigned int column = 0; column < numColumnsSubField; column++)
			{
				EXPECT_NEAR(expected(row, column), detector.test_getSubfield()(row, column), 0.015 * std::max(expected(row, column), detector.test_getSubfield()(row, column)));
			}
		}

		// Bias register map: check dimensions and content (unaltered)

		ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
		ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

		EXPECT_TRUE(arma::all(arma::vectorise(biasMap) == arma::vectorise(detector.test_getBiasRegisterMap())));

		// Smearing map: check dimensions and content (unaltered)

		ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
		ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

		EXPECT_TRUE(arma::all(arma::vectorise(smearingMap) == arma::vectorise(detector.test_getSmearingMap())));
	}
}











TEST_F(DetectorTest, applyOpenShutterSmearing)
{
	LOG_STARTING_OF_TEST


	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const double readoutTime = configParams.getDouble("CCD/ReadoutTime");	// 2
	const double exposureTime = readoutTime * 20.0;

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	// Initialise sub-pixel map, pixel map, bias register map, and smearing map

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::fmat(numRowsSubField, numColumnsSubField);
	subField.zeros();
	subField(0, 2) = 100.0;
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);



	// Open-shutter smearing

	detector.test_applyOpenShutterSmearing(exposureTime);



	// Sub-pixel map: check dimensions and content (unaltered)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap) == arma::vectorise(detector.test_getSubPixelMap())));

	// Pixel map: check dimensions and content (added open-shutter smearing)

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

	double expectedNoise = subField(0, 2) / numRowsSubField * readoutTime / exposureTime;

	EXPECT_TRUE(expectedNoise != 0.0);

	for(unsigned int row = 0; row < numRowsSubField; row++)
	{
		for(unsigned int column = 0; column < numColumnsSubField; column++)
		{
			if(column != 2)
			{
				EXPECT_EQ(0.0, detector.test_getSubfield()(row, column));
			}

			else
			{
				EXPECT_EQ(subField(row, column) + expectedNoise, detector.test_getSubfield()(row, column));
			}
		}
	}

	// Bias register map: check dimensions and content (added normal distribution)

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(biasMap) == arma::vectorise(detector.test_getBiasRegisterMap())));

	// Smearing map: check dimensions and content (unaltered)

	EXPECT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	EXPECT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

	for(unsigned int row = 0; row < numBiasPreScanRows; row++)
	{
		for(unsigned int column = 0; column < numColumnsSubField; column++)
		{
			if(column != 2)
			{
				EXPECT_FLOAT_EQ(smearingMap(row, column), detector.test_getSmearingMap()(row, column));
			}

			else
			{
				EXPECT_FLOAT_EQ(smearingMap(row, column) + expectedNoise, detector.test_getSmearingMap()(row, column));
			}
		}
	}
}









/**
 * Readout noise.
 *
 * Readout noise must be added to the pixel map and the bias register map.
 */
TEST_F(DetectorTest, addReadoutNoise)
{
	LOG_STARTING_OF_TEST

	// Construction

	configParams.setParameter("SubField/NumRows", "100");
	configParams.setParameter("SubField/NumColumns", "100");

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	const int readoutNoise = configParams.getInteger("CCD/ReadoutNoise");


	// Initialise sub-pixel map, pixel map, bias register map, and smearing map

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
	subField *= 1000.0;
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);



	// Readout noise

	detector.test_addReadoutNoise();



	// Readout noise distribution

	EXPECT_EQ(0.0, detector.test_getReadoutNoiseDistribution().mean());
	EXPECT_EQ(readoutNoise, detector.test_getReadoutNoiseDistribution().stddev());

	// Sub-pixel map: check dimensions and content (unaltered)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap) == arma::vectorise(detector.test_getSubPixelMap())));

	// Pixel map: check dimensions and content (added normal distribution)

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);
	arma::fmat residualSubField = detector.test_getSubfield() - subField;
	double stdDev = sqrt(arma::accu(residualSubField % residualSubField) / (numRowsSubField * numColumnsSubField));

	EXPECT_NEAR(0.0, mean(mean(residualSubField)), 0.6);
	EXPECT_NEAR(readoutNoise, stdDev, 0.01 * readoutNoise);

	// Bias register map: check dimensions and content (added normal distribution)

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	arma::fmat residualBiasRegisterMap = detector.test_getBiasRegisterMap() - biasMap;
	stdDev = sqrt(arma::accu(residualBiasRegisterMap % residualBiasRegisterMap) / (numBiasPreScanRows * numColumnsSubField));

	ASSERT_NEAR(0.0, mean(mean(residualBiasRegisterMap)), 0.6);
	ASSERT_NEAR(readoutNoise, stdDev, 0.05 * readoutNoise);

	// Smearing map: check dimensions and content (unaltered)

	EXPECT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	EXPECT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(smearingMap) == arma::vectorise(detector.test_getSmearingMap())));
}









/**
 * Gain.
 *
 * The values in the pixel map, bias register map, and smearing map must be divided by the
 * gain, on order to convert from [e- / pixel] to [ADU / pixel].
 */
TEST_F(DetectorTest, applyGain)
{
	LOG_STARTING_OF_TEST

	// Construction

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	const int gain = configParams.getInteger("CCD/Gain");

	// Initialise sub-pixel map, pixel map, bias register map, and smearing map

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);



	// Gain

	detector.test_applyGain();



	// Sub-pixel map: check dimensions and content (unaltered)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap) == arma::vectorise(detector.test_getSubPixelMap())));

	// Pixel map: check dimensions and content (divided by gain)

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subField / gain) == arma::vectorise(detector.test_getSubfield())));

	// Bias register map: check dimensions and content (divided by gain)

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(biasMap / gain) == arma::vectorise(detector.test_getBiasRegisterMap())));

	// Smearing map: check dimensions and content (divided by gain)

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(smearingMap / gain) == arma::vectorise(detector.test_getSmearingMap())));
}









/**
 * Electronic offset.
 *
 * The electronic offset must be added to the pixel map, bias register map, and smearing map.
 */
TEST_F(DetectorTest, addElectronicOffset)
{
	LOG_STARTING_OF_TEST

	// Construction

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	const int electronicOffset = configParams.getInteger("CCD/ElectronicOffset");

	// Initialise sub-pixel map, pixel map, bias register map, and smearing map

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);
	detector.test_setSubPixelMap(subPixelMap);

	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);
	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);
	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);
	detector.test_setSmearingMap(smearingMap);



	// Electronic offset

	detector.test_addElectronicOffset();



	// Sub-pixel map: check dimensions and content (unaltered)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap) == arma::vectorise(detector.test_getSubPixelMap())));

	// Pixel map: check dimensions and content (added electronic offset)

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subField + electronicOffset) == arma::vectorise(detector.test_getSubfield())));

	// Bias register map: check dimensions and content (added electronic offset)

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(biasMap + electronicOffset) == arma::vectorise(detector.test_getBiasRegisterMap())));

	// Smearing map: check dimensions and content (added electronic offset)

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(smearingMap + electronicOffset) == arma::vectorise(detector.test_getSmearingMap())));
}









/**
 * Digital saturation.
 */
TEST_F(DetectorTest, applyDigitalSaturation)
{
	LOG_STARTING_OF_TEST

	// Construction

	JitterFromRedNoise jitterGenerator(configParams);
	Platform platform(configParams, hdf5File, jitterGenerator);
	Sky sky(configParams);
	Telescope telescope(configParams, hdf5File, platform);
	Camera camera(configParams, hdf5File, telescope, sky);
	MyDetector detector(configParams, hdf5File, camera);

	// Configuration parameters

	const int numRowsSubField = configParams.getInteger("SubField/NumRows");
	const int numColumnsSubField = configParams.getInteger("SubField/NumColumns");

	const int numBiasPreScanRows = configParams.getInteger("SubField/NumBiasPrescanRows");
	const int numSmearingOverScanRows = configParams.getInteger("SubField/NumSmearingOverscanRows");

	const int numSubPixels = configParams.getInteger("SubField/SubPixels");

	const int digitalSaturationLimit = configParams.getInteger("CCD/DigitalSaturation");

	// Initialise the sub-pixel map, pixel map, bias register map, and smearing map
	// (make sure some pixel values exceed the digital saturation limit)

	int row;
	int column;

	std::default_random_engine engine;

	arma::fmat subPixelMap = arma::randu<arma::fmat>(numRowsSubField * numSubPixels, numColumnsSubField * numSubPixels);

	for(unsigned int index = 0; index < 5; index++)
	{
		std::uniform_real_distribution<double> randomRow(0, numRowsSubField * numSubPixels - 1);
		std::uniform_real_distribution<double> randomColumn(0, numColumnsSubField * numSubPixels - 1);

		row = randomRow(engine);
		column = randomColumn(engine);

		subPixelMap(row, column) = digitalSaturationLimit + 1;
	}

	detector.test_setSubPixelMap(subPixelMap);


	arma::fmat subField = arma::randu<arma::fmat>(numRowsSubField, numColumnsSubField);

	for(unsigned int index = 0; index < 5; index++)
	{
		std::uniform_real_distribution<double> randomRow(0, numRowsSubField - 1);
		std::uniform_real_distribution<double> randomColumn(0, numColumnsSubField - 1);

		row = randomRow(engine);
		column = randomColumn(engine);

		subField(row, column) = digitalSaturationLimit + 1;
	}

	detector.test_setSubfield(subField);

	arma::fmat biasMap = arma::randu<arma::fmat>(numBiasPreScanRows, numColumnsSubField);

	for(unsigned int index = 0; index < 5; index++)
	{
		std::uniform_real_distribution<double> randomRow(0, numBiasPreScanRows - 1);
		std::uniform_real_distribution<double> randomColumn(0, numColumnsSubField - 1);

		row = randomRow(engine);
		column = randomColumn(engine);

		biasMap(row, column) = digitalSaturationLimit + 1;

	}

	detector.test_setBiasRegisterMap(biasMap);

	arma::fmat smearingMap = arma::randu<arma::fmat>(numSmearingOverScanRows, numColumnsSubField);

	for(unsigned int index = 0; index < 5; index++)
	{
		std::uniform_real_distribution<double> randomRow(0, numSmearingOverScanRows - 1);
		std::uniform_real_distribution<double> randomColumn(0, numColumnsSubField - 1);

		row = randomRow(engine);
		column = randomColumn(engine);

		smearingMap(row, column) = digitalSaturationLimit + 1;

	}

	detector.test_setSmearingMap(smearingMap);



	// Digital saturation

	detector.test_applyDigitalSaturation();



	// Sub-pixel map: check dimensions and content (unaltered)

	ASSERT_EQ(numRowsSubField * numSubPixels, detector.test_getSubPixelMap().n_rows);
	ASSERT_EQ(numColumnsSubField * numSubPixels, detector.test_getSubPixelMap().n_cols);

	EXPECT_TRUE(arma::all(arma::vectorise(subPixelMap) == arma::vectorise(detector.test_getSubPixelMap())));

	// Pixel map: check dimensions and content (topped off at digital saturation limit)

	ASSERT_EQ(numRowsSubField, detector.test_getSubfield().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSubfield().n_cols);

	ASSERT_EQ(digitalSaturationLimit, detector.test_getSubfield().max());

	for(unsigned int row = 0; row < numRowsSubField; row++)
	{
		for(unsigned int column = 0; column < numColumnsSubField; column++)
		{
			double expected = std::min((float) digitalSaturationLimit, subField.at(row, column));
			ASSERT_EQ(expected, detector.test_getSubfield()(row, column));
		}
	}


	// Bias register map: check dimensions and content (topped off at digital saturation limit)

	ASSERT_EQ(numBiasPreScanRows, detector.test_getBiasRegisterMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getBiasRegisterMap().n_cols);

	ASSERT_EQ(digitalSaturationLimit, detector.test_getBiasRegisterMap().max());

	for(unsigned int row = 0; row < numBiasPreScanRows; row++)
	{
		for(unsigned int column = 0; column < numColumnsSubField; column++)
		{
			double expected = std::min((float) digitalSaturationLimit, biasMap.at(row, column));
			ASSERT_EQ(expected, detector.test_getBiasRegisterMap()(row, column));
		}
	}

	// Smearing map: check dimensions and content (topped off at digital saturation limit)

	ASSERT_EQ(numSmearingOverScanRows, detector.test_getSmearingMap().n_rows);
	ASSERT_EQ(numColumnsSubField, detector.test_getSmearingMap().n_cols);

	ASSERT_EQ(digitalSaturationLimit, detector.test_getSmearingMap().max());

	for(unsigned int row = 0; row < numSmearingOverScanRows; row++)
	{
		for(unsigned int column = 0; column < numColumnsSubField; column++)
		{
			double expected = std::min((float) digitalSaturationLimit, smearingMap.at(row, column));
			ASSERT_EQ(expected, detector.test_getSmearingMap()(row, column));
		}
	}
}









TEST_F(DetectorTest, convolveWithPsf)
{
	LOG_STARTING_OF_TEST

}









TEST_F(DetectorTest, getPlanarFocalPlaneCoordinatesOfSubfieldCorners)
{
	LOG_STARTING_OF_TEST

}









TEST_F(DetectorTest, getSolidAngleOfOnePixel)
{
	LOG_STARTING_OF_TEST
}
