#include "detector.h"

/**
 * Constructor.  Creates a detector object, based on the given configuration
 * parameters and attaches it to the given camera.
 *
 * @param configurationParameters: Configuration parameters for the detector.
 * @param camera:                  Camera to which to attach the detector.
 */

Detector::Detector()
{

	// Parse the parameters from the configuration file

	// Associate the camera

	// Allocate memory for the sub-pixel map, the pixel map, the bias register map,
	// and the smearing map

	initSubPixelMap();
	initPixelMap();
	initBiasMap();
	initSmearingMap();

	// Initialise the flatfield map and the CTE map

	initFlatfieldMap();
	initCteMap();

	// Random number generators

	photonNoiseGenerator.seed(photonNoiseSeed);
	readoutNoiseGenerator.seed(readoutNoiseSeed);
}








/**
 * Destructor.
 *
 * @post De-allocated memory of sub-field.
 * @post De-allocated memory of CTE map.
 * @post De-allocated memory of flatfield map.
 * @post De-allocated memory of bias register map.
 * @post De-allocated memory of smearing map.
 */

Detector::~Detector()
{

	// De-allocate the sub-pixel map

	for (unsigned int row = 0; row < numRowsSubPixelMap; row++)
	{
		delete[] subPixelMap[row];
	}

	delete[] subPixelMap;

	// De-allocate the flatfield map

	for (unsigned int row = 0; row < numRowsSubField * numSubPixelsPerPixel;
			row++)
	{
		delete[] flatfieldMap[row];
	}

	delete[] flatfieldMap;

	// De-allocate the pixel map and the CTE map

	for (unsigned int row = 0; row < numRowsSubField; row++)
	{
		delete[] pixelMap[row];
		delete[] cteMap[row];
	}

	delete[] pixelMap;
	delete[] cteMap;

	// De-allocate the bias register map

	for (unsigned int row = 0; row < numRowsBiasMap; row++)
	{
		delete[] biasMap[row];
	}

	delete[] biasMap;

	// De-allocate the smearing map

	for (unsigned int row = 0; row < numRowsSmearingMap; row++)
	{
		delete[] smearingMap[row];
	}

	delete[] smearingMap;
}





/**
 * Method that allocates memory for the sub-pixel map.
 */
void Detector::initSubPixelMap()
{
	// Allocate memory for the sub-pixel map

	subPixelMap = new double*[numRowsSubPixelMap];

	for (unsigned int row = 0; row < numRowsSubPixelMap; row++)
	{
		subPixelMap[row] = new double[numColumnsSubPixelMap];
	}
}








/**
 * Method that allocates memory for the pixel map.
 */
void Detector::initPixelMap()
{
	// Allocate memory for the pixel map

	pixelMap = new double*[numRowsSubField];

	for (unsigned int row = 0; row < numRowsSubField; row++)
	{
		pixelMap[row] = new double[numColumnsSubField];
	}
}







/**
 * Method that allocates memory for the bias register map.
 */
void Detector::initBiasMap()
{
	// Allocate memory for the bias register map

	biasMap = new double*[numRowsBiasMap];

	for (unsigned int row = 0; row < numRowsBiasMap; row++)
	{
		biasMap[row] = new double[numColumnsSubField];
	}
}






/**
 * Method that allocates memory for the smearing map.
 */
void Detector::initSmearingMap()
{
	// Allocate memory for the smearing map

	smearingMap = new double*[numRowsSmearingMap];

	for(unsigned int row = 0; row < numRowsSmearingMap; row++)
	{
		smearingMap[row] = new double[numColumnsSubField];
	}
}







/**
 * Method that allocates memory for the flatfield map and initialises it.
 */
void Detector::initFlatfieldMap()
{
	// Random number generation

	mt19937 flatfieldGenerator(flatfieldSeed);
	normal_distribution<double> flatfieldDistribution(0.0, 1.0);

	// Create a square map, filled with zeroes, in which the whole sub-field fits
	// at pixel level and for which the dimensions are a power of 2

	unsigned int dimensionPowerOf2 = 2;
	unsigned int maxSubFieldDimension = max(numRowsSubField,
			numColumnsSubField);

	while (dimensionPowerOf2 <= maxSubFieldDimension)
	{
		dimensionPowerOf2 *= 2;
	}

	double** flatFieldMapPowerOf2Dimensions = new double*[dimensionPowerOf2];

	for (unsigned int row = 0; row < dimensionPowerOf2; row++)
	{
		flatFieldMapPowerOf2Dimensions[row] = new double[dimensionPowerOf2];

		for (unsigned int column = 0; column < dimensionPowerOf2; column++)
		{
			flatFieldMapPowerOf2Dimensions[row][column] = 0.0;
		}
	}

	// Add variations at all spatial frequencies
	// Recursive process (base case: n = dimension of the map / 2)
	//		- add the same (random) value to pixels in the same (n,n) block
	//		- n /= 2
	//		- continue as long as n > 2

	unsigned int numBlocks = 2;
	double variation;

	// Loop over all block sizes: dimensionPowerOf2 / 2, dimensionPowerOf2 / 4,
	// dimensionPowerOf2 / 8,..., 2

	for (unsigned int blockSize = dimensionPowerOf2 / 2; blockSize >= 2;
			blockSize /= 2)
	{
		// Loop over all blocks of the current size

		for (unsigned int blockRow = 0; blockRow < numBlocks; blockRow++)
		{
			for (unsigned int blockColumn = 0; blockColumn < numBlocks;
					blockColumn++)
			{
				variation = flatfieldDistribution(flatfieldGenerator);

				for (unsigned row = blockRow * blockSize;
						blockRow < (blockRow + 1) * blockSize; row++)
				{
					for (unsigned column = blockColumn * blockSize;
							column < (blockColumn + 1) * blockSize; column++)
					{
						flatFieldMapPowerOf2Dimensions[row][column] +=
								variation;
					}
				}
			}
		}

		numBlocks *= 2;
	}

	// Normalise and subtract 0.5 -> [-0.5, 0.5]
	// Multiply by peak-to-peak noise amplitude

	double minValue = std::numeric_limits<double>::max();
	double maxValue = -std::numeric_limits<double>::max();

	for (unsigned int row = 0; row < dimensionPowerOf2; row++)
	{
		for (unsigned int column = 0; column < dimensionPowerOf2; column++)
		{
			minValue = std::min(minValue,
					flatFieldMapPowerOf2Dimensions[row][column]);// Look for the minimum value
			maxValue = std::max(maxValue,
					flatFieldMapPowerOf2Dimensions[row][column]);// Look for the maximum value
		}
	}

	for (unsigned int row = 0; row < dimensionPowerOf2; row++)
	{
		for (unsigned int column = 0; column < dimensionPowerOf2; column++)
		{
			flatFieldMapPowerOf2Dimensions[row][column] =
					((flatFieldMapPowerOf2Dimensions[row][column] - minValue)
							/ maxValue - 0.5)
							/ flatfieldPeak2PeakNoiseAmplitude;
		}
	}

	// Allocate memory for the flatfield map (at sub-pixel level)

	unsigned int numRowsFlatfieldMap = numSubPixelsPerPixel * numRowsSubField;
	unsigned int numColumnsFlatfieldMap = numSubPixelsPerPixel
			* numColumnsSubField;

	flatfieldMap = new double*[numRowsFlatfieldMap];

	for (unsigned int row = 0; row < numRowsFlatfieldMap; row++)
	{
		flatfieldMap[row] = new double[numColumnsFlatfieldMap];

		for (unsigned int column = 0; column < numColumnsFlatfieldMap; column++)
		{
			flatfieldMap[row][column] = 0.0;
		}
	}

	// The central part of the pixels has a sensitivity loss of less than 5%

	unsigned int edge = (int) ceil(
			numSubPixelsPerPixel * flatfieldIntraPixelWidth / 100.);

	int flatfieldRow, flatfieldColumn;

	// Loop over all pixels

	for (unsigned int pixelRow = 0; pixelRow < numRowsSubField; pixelRow++)
	{
		for (unsigned int pixelColumn = 0; pixelColumn < numColumnsSubField;
				pixelColumn++)
		{
			// Loop over all sub-pixels in the current pixel

			for (unsigned int row = 0; row < numSubPixelsPerPixel; row++)
			{
				for (unsigned int column = 0; column < numSubPixelsPerPixel;
						column++)
				{
					flatfieldRow = (pixelRow * numSubPixelsPerPixel) + row;
					flatfieldColumn = (pixelColumn * numSubPixelsPerPixel)
							+ column;

					// Edge: sensitivity loss of 5%

					if ((row < edge) || (column < edge)
							|| (column >= numSubPixelsPerPixel - edge)
							|| (row >= numSubPixelsPerPixel - edge))
					{
						flatfieldMap[flatfieldRow][flatfieldColumn] =
								(1
										+ flatFieldMapPowerOf2Dimensions[pixelRow][pixelColumn]
										+ flatfieldDistribution(
												flatfieldGenerator)
												* flatfieldWhiteNoise)
										* intraPixelSensitivity;
					}

					// Central part: no sensitivity loss

					else
					{
						flatfieldMap[flatfieldRow][flatfieldColumn] =
								(1
										+ flatFieldMapPowerOf2Dimensions[pixelRow][pixelColumn]
										+ flatfieldDistribution(
												flatfieldGenerator)
												* flatfieldWhiteNoise);
					}
				}
			}
		}
	}

	// De-allocate memory

	for(unsigned int row = 0; row < dimensionPowerOf2; row++)
	{
		delete[] flatFieldMapPowerOf2Dimensions[row];
	}

	delete[] flatFieldMapPowerOf2Dimensions;
}






/**
 * Method that allocates memory for the CTE map and initialises it.
 */
void Detector::initCteMap()
{
//	discrete_distribution<int> cteRowDistribution(0.0,
//			(double) numRowsSubField);
//	discrete_distribution<int> cteColumnDistribution(0.0,
//			(double) numColumnsSubField);
//
//	mt19937 cteRowGenerator(cteMapSeedRow);
//	mt19937 cteColumnGenerator(cteMapSeedColumn);

	for (unsigned int row = 0; row < numRowsSubField; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			cteMap[row][column] = meanCte;
		}
	}
}








/**
 * Method that resets the sub-field, the bias register map, and the smearing map
 * (i.e. all values are set to zero).
 *
 * @pre Sub-pixel map possibly filled with values from previous exposure.
 * @pre Pixel map possibly filled with values from previous exposure.
 * @pre Bias register map possibly filled with values from previous exposure.
 * @pre Smearing map possibly filled with values from previous exposure.
 *
 * @post Sub-pixel map filled with zeroes.
 * @post Pixel map filled with zeroes
 * @post Bias register map  filled with zeroes.
 * @post Smearing map filled with zeroes.
 */

void Detector::reset()
{

	// Set elements of sub-pixel map to zero

	for (unsigned int row = 0; row < numRowsSubPixelMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubPixelMap; column++)
		{
			subPixelMap[row][column] = 0.0;
		}
	}


	// Set elements of pixel map to zero

	for (unsigned int row = 0; row < numRows; row++)
	{
		for (unsigned int column = 0; column < numColumns; column++)
		{
			pixelMap[row][column] = 0.0;
		}
	}


	// Set elements of bias map to zero

	for (unsigned int row = 0; row < numRowsBiasMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			biasMap[row][column] = 0.0;
		}
	}


	// Set elements of smearing map to zero

	for (unsigned int row = 0; row < numRowsSmearingMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			smearingMap[row][column] = 0.0;
		}
	}
}









/**
 * Method that takes an exposure with the detector starting at the given time.
 *
 * The light is integrated during the given exposure time, during which the
 * detector experiences the effects of jitter and telescope drift (i.e. thermo-
 * elastic variations).  The background is assumed uniform for the whole
 * sub-field.
 *
 * Afterwards, the collected light is read out, convolving the image with the
 * PSF of the camera and adding various noise effects.
 *
 * @param startTime: Starting time of the exposure [s].
 * @param exposureTime: Duration of the exposure [s].
 *
 * @pre Sub-pixel map possibly filled with values from previous exposure.
 * @pre Pixel map possibly filled with values from previous exposure.
 * @pre Bias register map possibly filled with values from previous exposure.
 * @pre Smearing map possibly filled with values from previous exposure.
 *
 * @post Pixel unit in the pixel map: [ADU]
 * @post Pixel unit in the bias register map: [ADU]
 * @post Pixel unit in the smearing map: [ADU]
 */

void Detector::takeExposure(double startTime, double exposureTime)
{

	// Integration of point sources and background, taking into account jitter + drift.

	integrateLight(startTime, exposureTime);

	// Readout: convolution with PSF + noise effects

	readOut(exposureTime);
}









/**
 * Method that applies jitter to the exposure.  The jitter position are either
 * read from an input file or can be calculated using the input parameters which
 * are read from the configuration file.
 *
 * For each of the jitter position in the exposure, each star position is
 * re-calculated and the flux is summed in the sub-pixel map, which introduces
 * blurring.  Note that the convolution with the PSF is not yet done here.
 *
 * This method also adds the contribution of the background (zodiacal + galactic)
 * to the sub-pixel map.
 *
 * @param startTime: Starting time of the exposure for which jitter must be applied [s].
 *
 * @pre Sub-pixel map possibly filled with values from previous exposure.
 * @pre Pixel map possibly filled with values from previous exposure.
 * @pre Bias register map possibly filled with values from previous exposure.
 * @pre Smearing map possibly filled with values from previous exposure.
 *
 * @post Pixel unit of the sub-pixel map: [photons]
 * @post Re-calculated star position for each jitter step and summed the flux
 *       in the sub-pixel map, which introduces blurring. Also the contribution
 *       of the background (zodiacal + galactic) is added to the sub-pixel map.
 * @post Pixel map filled with zeroes.
 * @post Bias register map filled with zeroes.
 * @post Smearing map filled with zeroes.
 */

void Detector::integrateLight(double startTime, double exposureTime)
{

	// Reset the sub-field (i.e. get rid of the previous exposure, by zeroing the entire sub-field)

	reset();

	// Integration (incl. jitter) + background

	//camera.exposeSubField(this, exposureTime);

	// Apply flatfield (at sub-pixel level)

	applyFlatfield();

	// Rebin

	rebin();
}










/**
 * Method that adds the given flux value to the value of the sub-pixel that
 * corresponds to the given coordinates in the focal plane.  The flux value has
 * already been multiplied with the transmission efficiency but not with the
 * qunatum efficiency.  Also the exposure time has been taken into account already.
 *
 * @param rowFocalPlane: Row coordinate of the sub-pixel in the focal plane [mm].
 * @param yCoords: Column coordinate of the sub-pixel in the focal plane [mm].
 * @param flux: Flux to add to the sub-pixel map [photons].
 *
 * @pre Flux received at the given position in the focal plane not added to the
 *      sub-pixel map yet.
 * @pre Pixel map filled with zeroes.
 * @pre Bias register map filled with zeroes.
 * @pre Smearing map filled with zeroes.
 *
 * @post Flux received at the given position in the focal plane added to the sub-pixel map.
 * @post Pixel unit in the sub-pixel map: [photons]
 * @post Pixel map filled with zeroes.
 * @post Bias register map filled with zeroes.
 * @post Smearing map filled with zeroes.
 */
void Detector::addFlux(double rowFocalPlane, double columnFocalPlane,
		double flux)
{

	// Detector origin offset (pixel level)

	double rowOffset = (rowFocalPlane - originOffsetRow) / pixelSize;
	double columnOffset = (columnFocalPlane - originOffsetColumn) / pixelSize;

	// Detector orientation (pixel level)

	double column = columnOffset * cos(orientationAngle)
			- rowOffset * sin(orientationAngle);
	double row = columnOffset * sin(orientationAngle)
			+ rowOffset * cos(orientationAngle);

	// Sub-field incl. edge pixels (also correct for sub-field zeropoint)

	column = (column - subFieldZeroPointColumn + numEdgePixels)
			* numSubPixelsPerPixel;
	row = (row - subFieldZeroPointRow + numEdgePixels) * numSubPixelsPerPixel;

	// Add flux in this->subPixelMap at (row, column)

	if (this->isInSubPixelMap(row, column))
	{
		subPixelMap[(int) round(row)][(int) round(column)] += flux;
	}
}








/**
 * Method that checks whether the given (row, column) coordinates are in the
 * sub-pixel map.
 *
 * @param row:    Row coordinate.
 * @param column: Column coordinate.
 *
 * @return True if the given (row, column) coordinates are in the sub-pixel map;
 *         false otherwise.
 */

bool Detector::isInSubPixelMap(double row, double column)
{
	return (column >= 0) && (row >= 0) && (column < numColumnsSubPixelMap)
			&& (row < numRowsSubPixelMap);
}










/**
 * Method that adds the given flux value to (all sub-pixels of) the sub-pixel map.
 *
 * @param flux: Flux to add to the sub-pixel map [photons], already multiplied
 *        with the transmission efficiency of the camera but not with the
 *        quantum efficiency of the detector.  Also the exposure time has been
 *        taken into account already.
 *
 * @pre Given flux value not added to the sub-pixel map yet.
 * @pre Pixel map filled with zeroes.
 * @pre Bias register map filled with zeroes.
 * @pre Smearing map filled with zeroes.
 *
 * @pre  Flux value added to the sub-pixel map.
 * @post Pixel unit in the sub-pixel map: [photons]
 * @post Pixel map filled with zeroes.
 * @post Bias register map filled with zeroes.
 * @post Smearing map filled with zeroes.
 */
void Detector::addFlux(double flux)
{
	for (unsigned int row = 0; row < numRowsSubPixelMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubPixelMap; column++)
		{
			subPixelMap[row][column] += flux;
		}
	}
}









/**
 * Method that multiplies the sub-pixel map with the flatfield map.  The central
 * part (i.e. all pixels except the edge pixels) are multiplied element-wise.
 *
 * @pre Pixel value in the sub-pixel map: [photons].
 * @pre Flatfield map at sub-pixel level, excl. edge pixels.
 * @pre Pixel map filled with zeroes.
 * @pre Bias register map filled with zeroes.
 * @pre Smearing map filled with zeroes.
 *
 * @post Central part (i.e. everything except the edge pixels) of the sub-pixel
 *       map is flatfielded.
 * @post Pixel value in the sub-pixel map: [photons].
 * @post Pixel map filled with zeroes.
 * @post Bias register map filled with zeroes.
 * @post Smearing map filled with zeroes.
 */
void Detector::applyFlatfield()
{
	unsigned int numEdgeSubPixels = numEdgePixels * numSubPixelsPerPixel;

	// Loop over all elements in the sub-pixel map, except the edge pixels

	for (unsigned int row = numEdgeSubPixels;
			row < numRowsSubPixelMap - numEdgeSubPixels; row++)
	{
		for (unsigned int column = numEdgeSubPixels;
				column < numColumnsSubPixelMap - numEdgeSubPixels; column++)
		{
			subPixelMap[row][column] *=
					flatfieldMap[row - numEdgeSubPixels][column
							- numEdgeSubPixels];
		}
	}
}





/**
 * Method that rebins the sub-pixel map to pixel level and crops the edge pixels
 * that were added on each side to account for the edge effect.
 *
 * @pre Pixel value in the sub-pixel map: [photons].
 * @pre Pixel map filled with zeroes.
 * @pre Bias register map filled with zeroes.
 * @pre Smearing map filled with zeroes.
 *
 * @post Pixel value in the sub-pixel map: [photons]
 * @post Pixel value in the sub-pixel map: [photons]
 * @post Bias register map filled with zeroes.
 * @post Smearing map filled with zeroes.
 */
void Detector::rebin()
{
	// Rebinning is simply done by adding all values of the subpixels per pixel.

	for (unsigned int row = 0; row < numRowsSubField; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			double sum = 0;

			for (unsigned int rowSubPixelMap = row * numSubPixelsPerPixel;
					rowSubPixelMap < (row + 1) * numSubPixelsPerPixel;
					rowSubPixelMap++)
			{
				for (unsigned int columnSubPixelMap = column
						* numSubPixelsPerPixel;
						columnSubPixelMap < (column + 1) * numSubPixelsPerPixel;
						columnSubPixelMap++)
				{
					sum += subPixelMap[rowSubPixelMap][columnSubPixelMap];
				}
			}

			pixelMap[row][column] = sum;
		}
	}
}






/**
 * Method that reads out the detector and applies the following effects:
 *   - quantum efficiency
 * 	 - photon noise
 *   - full-well saturation (i.e. blooming)
 * 	 - CTE
 *   - open-shutter smearing
 * 	 - readout noise
 * 	 - gain
 * 	 - electronic offset (i.e. bias)
 * 	 - digital saturation
 *
 * @param exposureTime: Exposure time [s].
 *
 * @pre Pixel unit in the pixel map: [photons]
 * @pre Bias register map filled with zeroes.
 * @pre Smearing map filled with zeroes.
 *
 * @post Pixel unit in the pixel map: [ADU]
 * @post Pixel unit in the bias register map: [ADU]
 * @post Pixel unit in the smearing map: [ADU]
 */
void Detector::readOut(double exposureTime)
{

	// Apply quantum efficiency
	// Pixel units before: [photons]
	// Pixel units after: [electrons]

	applyQuantumEfficiency();

	// Apply poisson distributed photon noise
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	if (doPhotonNoise)
	{
		addPhotonNoise();
	}

	// Apply full-well saturation. A pixel has a maximum capacity of electrons (the full well capacity).
	// If photons free more electrons, the pixel saturates, and the electrons flow in the pixels above and below in
	// the same column (potential barriers are smallest in that direction).
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	applyFullWellSaturation();

	// Simulate the effects of the Charge Transfer Inefficiency (CTI). When the
	// CCD is read out, row after row, a part of the charge is always left behind
	// which then dribbles into the trailing pixels. This causes each star to have
	// a small "tail". Only visible when the CTI = 1 - CTE is poor.
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	applyCte();

	// Apply the effects of readout smearing due to an open shutter. Because there is no shutter,
	// the pixels are still receiving photons from the sky, while they are being transfered towards
	// the readout register.

	applyOpenShutterSmearing();

	// Each time the amplifier reads out a pixel, a tiny bit of noise is added.
	// Add the readout noise.
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	addReadoutNoise();

	// Apply the gain, to increase the dynamic range of the detector.
	// Pixel units before: [electrons]
	// Pixel units after: [ADU]

	applyGain();

	// Take into account the bias level (i.e. add the constant "zero" level
	// introduced by the amplifier).
	// Pixel units before: [ADU]
	// Pixel units after: [ADU]

	addElectronicOffset();

	// Take into account digital saturation. If even after dividing by the gain
	// the number of ADUs in a pixel is still higher than the analogue-digital
	// converter (ADC) can represent with its fixed amount of bits, clip all
	// values that are too high to the saturation level of the ADC.
	// Pixel units before: [ADU]
	// Pixel units after: [ADU]

	applyDigitalSaturation();
}








/**
 * Method that applies that quantum efficiency to the pixel.  The pixel values
 * are multiplied by the quantum efficiency of the detector.
 *
 * @pre Pixel unit in the pixel map: [photons]
 * @pre Bias register map filled with zeroes.
 * @pre Smearing map filled with zeroes.
 *
 * @pre Pixel unit in the pixel map: [electrons]
 * @post Bias register map filled with zeroes.
 * @post Smearing map filled with zeroes.
 */

void Detector::applyQuantumEfficiency()
{
	for (unsigned int row = 0; row < numRowsSubField; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			pixelMap[row][column] *= quantumEfficiency;
		}
	}
}







/**
 * Method that adds photon noise (i.e. shot noise) to the pixel map.  This type
 * of noise occurs because of the discrete/quantised nature of the electric
 * charge carried by the electrons (when counting them as representatives of
 * photons hitting the detectors).  It follows a Poisson distribution and each
 * pixel is treated independently of the other pixels.
 *
 * @pre Pixel unit in the pixel map: [electrons]
 * @pre No smearing map
 * @pre No bias register map
 *
 * @post Pixel unit in the pixel map: [electrons]
 * @post Added photon noise to the pixel map
 * @post Pixel unit in the smearing map: [electrons]
 * @post Added photon noise to the smearing map
 * @post No bias register map
 */
void Detector::addPhotonNoise()
{
	for (unsigned int row = 0; row < numRowsSubField; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			photonNoiseDistribution = poisson_distribution<int>(pixelMap[row][column]);
			pixelMap[row][column] = photonNoiseDistribution(photonNoiseGenerator);
		}
	}

	// Add photon noise to the smearing map

	for (unsigned int row = 0; row < numRowsSmearingMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			photonNoiseDistribution = poisson_distribution<int>(smearingMap[row][column]);
			smearingMap[row][column] = photonNoiseDistribution(photonNoiseGenerator);
		}
	}
}










/**
 * Method that applies the effect of full-well saturation (i.e. blooming) to the
 * pixel map.  If a pixel receives more electrons than the full-well saturation
 * limit (expressed in [electrons / pixel]), the additional electrons flow evenly
 * distributed in positive and negative charge-transfer direction.  Electrons
 * reaching the edge of the CCD will not be detected.
 *
 * @pre Pixel unit in the pixel map: [electrons]
 * @pre Pixel unit in the smearing map: [electrons]
 * @pre No bias register map
 * @pre Full-well saturation limit expressed in [electrons]
 *
 * @post Pixel unit in the pixel map: [electrons]
 * @post Effect of full-well saturation (i.e. blooming) applied to the pixel map
 * @post Pixel unit in the smearing map: [electrons]
 * @post No bias register map
 */
void Detector::applyFullWellSaturation()
{
	double pixelValue, numExcessElectrons;

	unsigned int jmod;// Row coordinate where excess electrons are transferred from and to

	for (unsigned int row = 0; row < numRowsSubField; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			pixelValue = pixelMap[row][column];

			// If the full-well saturation limit has been exceeded, distribute
			// the electrons evenly in the wells above and below until the
			// saturation has disappeared (stay in the same column!)

			if (pixelValue > fullWellSaturationLimit)
			{
				// Transfer excess electrons down

				jmod = row;
				numExcessElectrons = (pixelValue - fullWellSaturationLimit)
						/ 2.0;    // Move half of the excess electrons down...

				while (numExcessElectrons > 0 && jmod < numRowsSubField)
				{
					pixelMap[jmod][column] -= numExcessElectrons;
					jmod++;

					// Electrons reaching the edge of the CCD will not be detected

					if (jmod < numRowsSubField)
					{
						pixelMap[jmod][column] += numExcessElectrons;

						// Make sure the pixel you move the excess electrons to
						// does not get saturated too

						if (pixelMap[jmod][column] > fullWellSaturationLimit)
						{
							numExcessElectrons = pixelMap[jmod][column]
									- fullWellSaturationLimit;
						}
					}
				}

				// Transfer excess electrons up

				jmod = row;
				numExcessElectrons = (pixelValue - fullWellSaturationLimit)
						/ 2.0;    // ...and the rest of the excess electrons up

				while (numExcessElectrons > 0 && jmod >= 0)
				{
					pixelMap[jmod][column] -= numExcessElectrons;
					jmod--;

					// Electrons reaching the edge of the CCD will not be detected

					if (jmod >= 0)
					{
						pixelMap[jmod][column] += numExcessElectrons;

						// Make sure the pixel you move the excess electrons to does not get saturated too

						if (pixelMap[jmod][column] > fullWellSaturationLimit)
						{
							numExcessElectrons = pixelMap[jmod][column]
									- fullWellSaturationLimit;
						}
					}
				}
			}
		}
	}
}









/**
 * Method that applies the effect of the charge-transfer (in)efficiency to the
 * pixel map. Assumed is that the serial register has a CTE of 1, unlike the
 * CCD that have a CTE map.
 *
 * @pre Pixel unit in the pixel map: [electrons]
 * @pre Pixel unit in the smearing map: [electrons]
 * @pre No bias register map
 *
 * @post Pixel unit in the pixel map: [electrons]
 * @post Applied the effect of the charge-transfer (in)efficiency of the CCD to
 *       the pixel map.
 * @post Pixel unit in the smearing map: [electrons]
 * @post No bias register map
 */
void Detector::applyCte()
{

}










/**
 * Method that applies the effect of readout smearing to the pixel map. This
 * effect is due to the fact that - because of the absence of a shutter (which
 * is common in space-based instruments) - the CCD still receives light during
 * frame transfer.  The flux of each pixel is affected by the flux of the pixels
 * in the same column.  Because the CCD is exposed during the whole readout and
 * multiple exposures are created, also the pixels further away from the readout
 * register have their influence.
 *
 * A smearing map is created and will be used in photometry to remove the
 * smearing effect from the pixel map.
 *
 * @pre Pixel unit in the pixel map: [electrons]
 * @pre Pixel unit in the smearing map: [electrons]
 * @pre No bias register map
 *
 * @post Pixel unit in the pixel map: [electrons]
 * @post Applied the effect of readout smearing to the pixel map
 * @post Pixel unit in the smearing map: [electrons]
 * @post Applied the effect of readout smearing to the smearing map
 * @post No bias register map
 */
void Detector::applyOpenShutterSmearing()
{

}









/**
 * Method that applies the readout noise to the pixel map and initialises the
 * bias register map.
 *
 * Readout noise occurs due to the imperfect nature of the CCD amplifiers.  When
 * the electrons are transferred to the amplifier, the induced voltage is
 * measured.  However, this measurement is not perfect, but gives a value which
 * is on average correct, with the readout noise as standard deviation.  So
 * readout noise is simply a measure of this scatter around the true value.
 * Its value is expressed in electrons because, after all, the packet of charge
 * is made up of electrons.
 *
 * @pre Pixel unit in the pixel map: [electrons]
 * @pre Pixel unit in the smearing map: [electrons]
 * @pre No bias register map
 * @pre Readout noise expressed in [electrons]
 *
 * @post Pixel unit in the pixel map: [electrons]
 * @post Added readout noise to the pixel map
 * @post Pixel unit in the smearing map: [electrons]
 * @post Pixel unit in the bias register map: [electrons]
 * @post Initialised the bias register map with readout noise
 */
void Detector::addReadoutNoise()
{
	readoutNoiseDistribution = normal_distribution<double>(0.0, readoutNoise);


	// Add readout noise to the pixel map

	for(unsigned int row = 0; row < numRowsSubField; row++)
	{
		for(unsigned int column = 0; column < numColumnsSubField; column++)
		{
			pixelMap[row][column] += readoutNoiseDistribution(readoutNoiseGenerator);
		}
	}

	// Initialise the bias  map with readout noise

	for(unsigned int row = 0; row < numRowsBiasMap; row++)
	{
		for(unsigned int column = 0; column < numColumnsSubField; column++)
		{
			biasMap[row][column] += readoutNoiseDistribution(readoutNoiseGenerator);
		}
	}
}








/**
 * Method that divides the bias map, smearing map, and pixel map by the
 * detector gain. This relates the number of electrons per pixels to the number
 * of counts (i.e. ADU) per pixel and converts these three maps from electrons
 * to ADU:
 *
 * @pre Pixel unit in the pixel map: [electrons]
 * @pre Pixel unit in the smearing map: [electrons]
 * @pre Pixel unit in the bias map: [electrons]
 *
 * @post Pixel unit in the pixel map: [ADU]
 * @post Pixel unit in the smearing map: [ADU]
 * @post Pixel unit in the bias  map: [ADU]
 */
void Detector::applyGain()
{
	// Multiply the pixel map with the gain

	for (unsigned int row = 0; row < numRowsSubField; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			pixelMap[row][column] *= gain;
		}
	}

	// Multiply the bias map with the gain

	for (unsigned int row = 0; row < numRowsBiasMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			biasMap[row][column] *= gain;
		}
	}

	// Multiply the smearing map with the gain

	for (unsigned int row = 0; row < numRowsSmearingMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			smearingMap[row][column] *= gain;
		}
	}
}








/**
 * Method that adds the electronic offset (i.e. bias level) to the pixel map,
 * smearing map, and bias register map to avoid negative readout values.
 *
 * @pre Pixel unit in the pixel map: [ADU]
 * @pre Pixel unit in the smearing map: [ADU]
 * @pre Pixel unit in the bias register map: [ADU]
 * @pre Electronic offset (i.e. bias level) expressed in [ADU]
 *
 * @post Pixel unit in the pixel map: [ADU]
 * @post Electronic offset added to the pixel map
 * @post Pixel unit in the smearing map: [ADU]
 * @post Electronic offset added to the smearing map
 * @post Pixel unit in the bias register map: [ADU]
 * @post Electronic offset added to the bias register map
 */
void Detector::addElectronicOffset()
{

	// Add the electronic offset to the pixel map

	for (unsigned int row = 0; row < numRowsSubField; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			pixelMap[row][column] += electronicOffset;
		}
	}

	// Add the electronic offset to the bias register map

	for (unsigned int row = 0; row < numRowsBiasMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			biasMap[row][column] += electronicOffset;
		}
	}

	// Add the electronic offset to the smearing map

	for (unsigned int row = 0; row < numRowsSmearingMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			smearingMap[row][column] += electronicOffset;
		}
	}
}








/**
 * Method that applies the effect of digital saturation to the pixel map,
 * smearing map, and bias register map.  This means that the pixel values in
 * these maps (expressed in [ADU / pixel]) are topped off to the digital saturation
 * limit of the detector (also expressed in [ADU / pixel]).
 *
 * @pre Pixel unit in the pixel map: [ADU]
 * @pre Pixel unit in the smearing map: [ADU]
 * @pre Pixel unit in the bias register map: [ADU]
 * @pre Digital saturation limit expressed in [ADU / pixel]
 *
 * @post Pixel unit in the pixel map: [ADU]
 * @post Values in the pixel map topped off at the digital saturation limit
 * @post Pixel unit in the smearing map: [ADU]
 * @post Values in the smearing map topped off at the digital saturation limit
 * @post Pixel unit in the bias register map: [ADU]
 * @post Values in the bias register map topped off at the digital saturation
 *       limit
 */
void Detector::applyDigitalSaturation()
{
	// Top off the values in the pixel map

	for (unsigned int row = 0; row < numRowsSubField; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			if (pixelMap[row][column] > digitalSaturationLimit)
			{
				pixelMap[row][column] = digitalSaturationLimit;
			}
		}
	}

	// Top off the values in the bias register map

	for (unsigned int row = 0; row < numRowsBiasMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			if (biasMap[row][column] > digitalSaturationLimit)
			{
				biasMap[row][column] = digitalSaturationLimit;
			}
		}
	}

	// Top off the values in the smearing map

	for (unsigned int row = 0; row < numRowsSmearingMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsSubField; column++)
		{
			if (smearingMap[row][column] > digitalSaturationLimit)
			{
				smearingMap[row][column] = digitalSaturationLimit;
			}
		}
	}
}
