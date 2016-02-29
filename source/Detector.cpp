#include "Detector.h"

#include "Units.h"

/**
 * Constructor.  
 *
 * \param hdf5file:                HFD5 file to write the detector images to.
 * \param configurationParameters: Configuration parameters for the detector.
 * \param camera:                  Camera to which to attach the detector.
 */
Detector::Detector(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera)
: HDF5Writer(hdf5file), internalTime(0.0), camera(camera), imageNr(0)
{
	// Create the groups in the HDF5 file where the different maps (i.e. pixel map,
	// bias register map, smearing map, etc.) will be saved. This needs to be done
	// BEFORE other methods write arrays to HDF5.

	initHDF5Groups();

	// Parse the parameters from the configuration file.

	configure(configParam);

	// Attach the camera

	// Allocate memory for the different maps

	pixelMap.zeros(numRowsPixelMap, numColumnsPixelMap);
	subPixelMap.zeros(numRowsSubPixelMap, numColumnsSubPixelMap);
	biasMap.zeros(numRowsBiasMap, numColumnsPixelMap);
	smearingMap.zeros(numRowsSmearingMap, numColumnsPixelMap);
	flatfieldMap.ones(numRowsSubPixelMap, numColumnsSubPixelMap);
//	cteMap.zeros(numRowsPixelMap, numColumnsPixelMap);

	// Generate the flatfield map 

	generateFlatfieldMap();

//	// Generate the CTE map
//
//	generateCteMap();

	// Set the seeds of the random number generators

	photonNoiseGenerator.seed(photonNoiseSeed);
	readoutNoiseGenerator.seed(readoutNoiseSeed);
}





/**
 * Destructor.
 *
 */
Detector::~Detector()
{

}











/**
 * \brief Configure the Detector object using the ConfigurationParameters
 * 
 * \param configParam: the configuration parameters 
 **/

 void Detector::configure(ConfigurationParameters &configParam)
 {
 	// Configuration parameters for the CCD detector

    originOffsetX           = configParam.getDouble("CCD/OriginOffsetX");
    originOffsetY           = configParam.getDouble("CCD/OriginOffsetY");
    orientationAngle        = deg2rad(configParam.getDouble("CCD/Orientation"));
    numRows                 = configParam.getInteger("CCD/NumRows");
    numColumns              = configParam.getInteger("CCD/NumColumns");
    pixelSize               = configParam.getDouble("CCD/PixelSize");
    gain                    = configParam.getInteger("CCD/Gain");
    quantumEfficiency       = configParam.getDouble("CCD/QuantumEfficiency");
    fullWellSaturationLimit = configParam.getLong("CCD/FullWellSaturation");
    digitalSaturationLimit  = configParam.getLong("CCD/DigitalSaturation");
    readoutNoise            = configParam.getDouble("CCD/ReadoutNoise");
    electronicOffset        = configParam.getInteger("CCD/ElectronicOffset");
    readoutTime             = configParam.getDouble("CCD/ReadoutTime");
    flatfieldNoiseAmplitude = configParam.getDouble("CCD/FlatfieldPtPNoise");
    meanCte                 = configParam.getDouble("CCD/CTEMean");

    // Configuration parameters for the subfield

    subFieldZeroPointRow    = configParam.getInteger("SubField/ZeroPointRow");
    subFieldZeroPointColumn = configParam.getInteger("SubField/ZeroPointColumn");
	numRowsPixelMap         = configParam.getInteger("SubField/NumRows");
	numColumnsPixelMap      = configParam.getInteger("SubField/NumColumns");
	numRowsBiasMap          = configParam.getInteger("SubField/NumBiasPrescanRows");
	numRowsSmearingMap      = configParam.getInteger("SubField/NumSmearingOverscanRows");
	numSubPixelsPerPixel    = configParam.getInteger("SubField/SubPixels");

    // Configuration parameters for the noise source random seeds

	readoutNoiseSeed        = configParam.getLong("RandomSeeds/ReadOutNoiseSeed");
	photonNoiseSeed         = configParam.getLong("RandomSeeds/PhotonNoiseSeed");
	flatfieldSeed           = configParam.getLong("RandomSeeds/FlatFieldSeed");
//	cteMapSeed              = configParam.getLong("RandomSeeds/CTESeed");

	// Derive the dimensions of the sub-pixel map

	numRowsSubPixelMap    = numRowsPixelMap    * numSubPixelsPerPixel;	// TODO Add edge pixels
	numColumnsSubPixelMap = numColumnsPixelMap * numSubPixelsPerPixel;	// TODO Add edge pixels

	numEdgePixels = 0;
 }











///**
// * \brief: Generate CTE map.  This map is generated at pixel level and currently
// *         the value of all elements in the CTE map are set to the mean CTE.
// *
// * NOTE: In a later version, we can introduce pixels and/or rows of pixels (in the
// *       pixel map) with a lower CTE, based on random distributions.
// */
//void Detector::generateCteMap()
//{
//	cteMap = meanCte;
//
//	// Random pixels with lower CTE
//
//	// Random rows of pixels with lower CTE
//}












/**
 * \brief: Generate the (random) flatfield variations.  This map is generated
 *		   at sub-pixel level but without the edge pixels.
 */
void Detector::generateFlatfieldMap()
{

	Log.info("Detector: generating flatfield map.");

	// Random number generation

	mt19937 flatfieldGenerator(flatfieldSeed);
	normal_distribution<double> flatfieldDistribution(0.0, 1.0);

	// Create a square map, filled with zeroes, in which the whole subfield fits,
	// and for which the dimensions are a power of 2

	unsigned int NrowsSquareMap = 2;
	unsigned int maxFlatfielMapDimension = max(flatfieldMap.n_rows, flatfieldMap.n_cols);

	while (NrowsSquareMap <= maxFlatfielMapDimension)
	{
		NrowsSquareMap *= 2;
	}

	arma::Mat<float> squareMap(NrowsSquareMap, NrowsSquareMap);
	squareMap.ones();

	// Add variations at all spatial frequencies
	// This is done by dividing the square map into:
	// 		overlapping blocks of size (N/2)x(N/2)
	// 		overlapping blocks of size (N/4)x(N/4)
	//      overlapping blocks of size (N/8)x(N/8)
	//      ...
	//
	// and add a gaussian noise to each of those blocks

	// Loop over all block sizes: N/2, N/4, N/8, ...

	for (unsigned int blockSize = NrowsSquareMap / 2; blockSize >= 2; blockSize /= 2)
	{
		// Loop over all overlapping blocks, and add a small variation

		for (unsigned int blockRow = 0; blockRow < NrowsSquareMap - blockSize; blockRow += blockSize)
		{
			for (unsigned int blockColumn = 0; blockColumn < NrowsSquareMap - blockSize; blockColumn += blockSize)
			{
				const double variation = flatfieldDistribution(flatfieldGenerator);
				squareMap(arma::span(blockRow, blockRow + blockSize - 1), arma::span(blockColumn, blockColumn + blockSize - 1)) += variation;
			}
		}
	}

	// Normalise and subtract 0.5 -> all values are in [-0.5, 0.5]
	// Multiply by peak-to-peak noise amplitude -> all values are in [-peakAmplitude/2, +peakAmplitude/2]

	double minValue = squareMap.min();
	double maxValue = squareMap.max();
	squareMap = ((squareMap - minValue) / (maxValue - minValue) - 0.5) * flatfieldNoiseAmplitude;

	// Copy a part of the squareMatrix corresponding to the size of the flatfieldMap, into the flatfieldMap

	flatfieldMap = squareMap.submat(0, 0, flatfieldMap.n_rows - 1, flatfieldMap.n_cols - 1);

	// Save the flatfield in the HDF5 file

	hdf5File.writeArray("/Flatfield", "flatfield", flatfieldMap);

}












/**
 * \brief: Zeroes the sub-pixel, pixel, bias register, and the smearing maps.
 *
 * \pre Sub-pixel, pixel, bias register, and smearing maps filled with values from previous exposure.
 *
 * \post Sub-pixel, pixel, bias register, and smearing maps filled with zeroes.
 */
void Detector::reset()
{
	subPixelMap.zeros();
	pixelMap.zeros();
	biasMap.zeros();
	smearingMap.zeros();
}









/**
 * \brief: Take an exposure with the detector starting at the given time.
 *		   The light is integrated during the given exposure time, during which 
 *         the detector experiences the effects of jitter and thermo-elastic telescope 
 *         drift. The background is assumed uniform for the whole subfield.
 *         Afterwards, the collected light is read out, convolving the image with the
 *         PSF of the camera and adding various noise effects.
 *
 * \param startTime:    Starting time of the exposure [s].
 * \param exposureTime: Duration of the exposure [s].
 * 
 * \return endTime:     Time after the exposure (startTime + exposureTime + readoutTime)
 *
 * \pre Sub-pixel, pixel, bias register, and smearing map filled with values from previous exposure.
 *
 * \post Pixel unit in the pixel, bias register, and smearing maps: [ADU]
 */

double Detector::takeExposure(double startTime, double exposureTime)
{
	// Advance the internal clock until the given start time

	internalTime = startTime;

	// Integration of point sources and background, taking into account jitter + drift.

	Log.debug("Detector: Integrating light for exposure " + to_string(imageNr) + " with exposure time = " + to_string(exposureTime));

	integrateLight(startTime, exposureTime);

	// Include noise effects like readout noise, photon noise, full well saturation, etc.
	// Note: readOut() needs the exposure time to compute the open shutter smearing.

	Log.debug("Detector: Adding noise effects to exposure " + to_string(imageNr));

	readOut(exposureTime);

	// Write the CCD subfield to the HDF5 file

	Log.debug("Detector: Writing PixelMap " + to_string(imageNr) + " to HDF5 file.");

	writePixelMapToHDF5();

	// Advance the internal clock

	internalTime += exposureTime + readoutTime;

	return internalTime;
}










/**
 * \brief: During an exposure, this method makes the detector integrates the light
 *         in small steps. During each step the slight change of star positions due
 *         to spacecraft jitter is taken into account. 
 *         Besides jitter, also the sky background, and the flatfield is taken into 
 *         account. The sub-pixel map is rebinned in a pixel map.
 *
 * NOTE: The convolution with the PSF is not yet done here.
 *
 * \param startTime: Starting time of the exposure for which jitter must be applied [s].
 *
 * \pre Sub-pixel, pixel, bias register, and smearing map filled with values from previous exposure.
 *
 * \post Pixel unit of the sub-pixel map: [photons].
 * \post Pixel, bias register, and smearing map filled with zeroes.
 */

void Detector::integrateLight(double startTime, double exposureTime)
{

	// Reset the sub-field (i.e. get rid of the previous exposure, by zeroing the entire sub-field)

	Log.debug("Detector: Resetting subfield array.");

	reset();

	// Integration (incl. jitter) + background

	camera.exposeDetector(*this, startTime, exposureTime);

	// Apply flatfield (at sub-pixel level)

	Log.debug("Detector: applying flatfield.");

	applyFlatfield();

	// Rebin from a subpixel map to a pixel map

	Log.debug("Detector: Rebinning subpixel map into pixel map.");

	rebin();
}












/**
 * \brief: Add the given flux value to the value of the sub-pixel that
 *         corresponds to the given coordinates in the focal plane.
 *
 * \param xFPprime   x-coordinate of the sub-pixel in the focal plane in the FP' reference frame [mm].
 * \param yFPprime   y-coordinate of the sub-pixel in the focal plane in the FP' reference frame [mm].
 * \param flux       Flux to add to the sub-pixel map [photons].
 *
 * \pre Pixel, bias register, and smearing maps filled with zeroes.
 *
 * \post Pixel unit in the sub-pixel map: [photons].
 * \post Pixel, bias register, and smearing maps filled with zeroes.
 */

void Detector::addFlux(double xFPprime, double yFPprime, double flux)
{

	// Detector origin offset (pixel level)

	double rowOffset    = (xFPprime - originOffsetY) / pixelSize;
	double columnOffset = (yFPprime - originOffsetX) / pixelSize;

	// Detector orientation (pixel level)

	double column = columnOffset * cos(orientationAngle) - rowOffset * sin(orientationAngle);
	double row    = columnOffset * sin(orientationAngle) + rowOffset * cos(orientationAngle);

	// Sub-field incl. edge pixels (also correct for sub-field zeropoint)

	column = (column - subFieldZeroPointColumn + numEdgePixels) * numSubPixelsPerPixel;
	row    = (row    - subFieldZeroPointRow    + numEdgePixels) * numSubPixelsPerPixel;

	// Add flux in this->subPixelMap at (row, column)

	if (isInSubPixelMap(row, column))
	{
		subPixelMap((int) round(row), (int) round(column)) += flux;
	}
}














/**
 * \brief Verify if a point with given planar focal plane coordinates is in the subfield
 * 
 * \param xFPprime    Planar focal plane x-coordinate in the FP' reference frame [mm]
 * \param yFPprime    Planar focal plane y-coordinate in the FP' reference frame [mm]
 * 
 * \return true if the point is in the subfield on the CCD, false otherwise.
 */

bool Detector::isInSubfield(const double xFPprime, const double yFPprime)
{
	// Convert to pixel coordinates in the unrotated CCD reference frame

	double rowUnrot = (xFPprime - originOffsetY) / pixelSize;
	double colUnrot = (yFPprime - originOffsetX) / pixelSize;

	// Compute the coordinates in the rotated CCD reference frame

	double colRot = colUnrot * cos(orientationAngle) - rowUnrot * sin(orientationAngle);
	double rowRot = colUnrot * sin(orientationAngle) + rowUnrot * cos(orientationAngle);

	// Check wether these pixel coordinates falls on the subfield

	return    (colRot >= subFieldZeroPointColumn) && (colRot < subFieldZeroPointColumn + numColumnsPixelMap)
	       && (rowRot >= subFieldZeroPointRow)    && (rowRot < subFieldZeroPointRow + numRowsPixelMap);
}













/**
 * \brief: Check whether the given (row, column) coordinates are in the
 *         sub-pixel map.
 *
 * NOTE: The input parameters row & column come from an coordinate transformation
 *       in the focal plane, and as a result are not necessarily integers. For the
 *       \brief of this function it's not necessary to round them to the nearest
 *       integer. 
 *
 * \param row:    Row coordinate     [sub-pixel].
 * \param column: Column coordinate  [sub-pixel].
 *
 * \return True if the given (row, column) coordinates are in the sub-pixel map;
 *         false otherwise.
 */
bool Detector::isInSubPixelMap(double row, double column)
{
	return (column >= 0) && (row >= 0) && (column < numColumnsSubPixelMap) && (row < numRowsSubPixelMap);
}










/**
 * \brief: Add the given flux value to (all sub-pixels of) the sub-pixel map.
 *
 * \param flux: Flux to add to the sub-pixel map [photons].
 *
 * \pre Pixel, bias register, and smearing maps filled with zeroes.
 *
 * \post Pixel unit in the sub-pixel map: [photons].
 * \post Pixel, bias register, and smearing maps filled with zeroes.
 */
void Detector::addFlux(double flux)
{
	subPixelMap += flux;
}










/**
 * \brief: Convolve the sub-pixel map with the given PSF, keeping the same dimensions.
 *
 * \param psf: PSF.
 */
void Detector::convolveWithPsf(arma::Mat<float> psf)
{
	subPixelMap = arma::conv2(subPixelMap, psf, "same");
}











/**
 * \brief: Multiply the sub-pixel map with the flatfield.
 * 
 * NOTE: The sub-pixel map contains extra edge pixels, but the flatfield
 *       map does not. These edge pixels are excluded from this flatfield
 *       multiplication.
 *
 * \pre Unit of the sub-pixels: [photons].
 * \pre Flatfield map at sub-pixel level, excl. edge pixels.
 * \pre Pixel, bias register, and smearing maps filled with zeroes.
 *
 * \post Pixel value in the sub-pixel map: [photons].
 * \post Pixel, bias, and smearing maps filled with zeroes.
 */
void Detector::applyFlatfield()
{
	const unsigned int numEdgeSubPixels = numEdgePixels * numSubPixelsPerPixel;
	const unsigned int beginRow = numEdgePixels;
	const unsigned int beginCol = numEdgePixels;
	const unsigned int endRow = numRowsSubPixelMap - numEdgeSubPixels - 1;
	const unsigned int endCol = numColumnsSubPixelMap - numEdgeSubPixels - 1;

	subPixelMap.submat(beginRow, beginCol, endRow, endCol) = subPixelMap.submat(beginRow, beginCol, endRow, endCol) % flatfieldMap;
}









/**
 * \brief: Rebin the sub-pixel map to pixel level and crop the edge pixels.
 *
 * \pre Unit of the pixel value in the sub-pixel map: [photons].
 * \pre Pixel, bias register, and smearing map filled with zeroes.
 *
 * \post Unit of pixel values in the sub-pixel map: [photons].
 * \post Bias register, and smearing maps filled with zeroes.
 */
void Detector::rebin()
{
	// Rebinning is simply done by adding all values of the sub-pixels per pixel.

	for (unsigned int row = 0; row < numRowsPixelMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsPixelMap; column++)
		{
			const unsigned int beginRow = row * numSubPixelsPerPixel;
			const unsigned int beginCol = column * numSubPixelsPerPixel;
			const unsigned int endRow = (row + 1) * numSubPixelsPerPixel - 1;
			const unsigned int endCol = (column + 1) * numSubPixelsPerPixel - 1;

			pixelMap(row, column) = arma::accu(subPixelMap.submat(beginRow, beginCol, endRow, endCol));
		}
	}
}












/**
 * \brief: Reads out the detector and apply the following effects:
 *   		- quantum efficiency
 * 	 		- photon noise
 *   		- full-well saturation (i.e. blooming)
 * 	 		- CTE
 *   		- open-shutter smearing
 * 	 		- readout noise
 * 	 		- gain
 * 	 		- electronic offset (i.e. bias)
 * 	 		- digital saturation
 *
 * \param exposureTime: Exposure time [s].
 *
 * \pre Pixel unit in the pixel map: [photons].
 * \pre Bias register and smearing maps filled with zeroes.
 *
 * \post Pixel unit in the pixel, bias register, and smearing maps: [ADU].
 */
void Detector::readOut(float exposureTime)
{

	// Apply quantum efficiency
	// Pixel units before: [photons]
	// Pixel units after: [electrons]

	applyQuantumEfficiency();

	// Apply poisson distributed photon noise
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	if (includePhotonNoise)
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

	applyOpenShutterSmearing(exposureTime);

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
 * \brief: Apply quantum efficiency to the pixel.  The pixel values
 *         are multiplied by the quantum efficiency of the detector.
 *
 * \pre Pixel unit in the pixel map: [photons].
 * \pre Bias and smearing maps filled with zeroes.
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \post Bias and smearing maps filled with zeroes.
 */
void Detector::applyQuantumEfficiency()
{
	Log.debug("Detector: applying quantum efficiency");

	pixelMap *= quantumEfficiency;
}










/**
 * \brief: Add photon noise (i.e. shot noise) to the pixel and smearing maps. 
 *         It follows a Poisson distribution and each pixel is treated 
 *         independently of the other pixels.
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre No bias register or smearing maps.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map.
 */
void Detector::addPhotonNoise()
{
	Log.debug("Detector: adding photon noise");

	// Add photon noise to the pixel map

	for (unsigned int row = 0; row < numRowsPixelMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsPixelMap; column++)
		{
			photonNoiseDistribution = poisson_distribution<long>(pixelMap(row, column));
			pixelMap(row, column) = photonNoiseDistribution(photonNoiseGenerator);
		}
	}

	// Add photon noise to the smearing map

	for (unsigned int row = 0; row < numRowsSmearingMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsPixelMap; column++)
		{
			photonNoiseDistribution = poisson_distribution<long>(smearingMap(row, column));
			smearingMap(row, column) = photonNoiseDistribution(photonNoiseGenerator);
		}
	}
}










/**
 * \brief: Apply the effect of full-well saturation (i.e. blooming) to the
 *         pixel map.  If a pixel receives more electrons than the full-well saturation
 *         limit (expressed in [electrons / pixel]), the additional electrons flow evenly
 *         distributed in positive and negative charge-transfer direction.  Electrons
 *         reaching the edge of the CCD will not be detected.
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre Pixel unit in the smearing map: [electrons].
 * \pre No bias register map.
 * \pre Full-well saturation limit expressed in [electrons].
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Effect of full-well saturation (i.e. blooming) applied to the pixel map.
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map.
 */
void Detector::applyFullWellSaturation()
{
	Log.debug("Detector: applying full well saturation");

	double pixelValue, numExcessElectrons;

	int jmod;// Row coordinate where excess electrons are transferred from and to

	for (unsigned int row = 0; row < numRowsPixelMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsPixelMap; column++)
		{
			pixelValue = pixelMap(row, column);

			// If the full-well saturation limit has been exceeded, distribute
			// the electrons evenly in the wells above and below until the
			// saturation has disappeared (stay in the same column!)

			if (pixelValue > fullWellSaturationLimit)
			{
				// Transfer excess electrons down

				jmod = row;
				numExcessElectrons = (pixelValue - fullWellSaturationLimit) / 2.0;   // Move half of the excess electrons down...

				while (numExcessElectrons > 0 && jmod < numRowsPixelMap)
				{
					pixelMap(jmod, column) -= numExcessElectrons;
					jmod++;

					// Electrons reaching the edge of the CCD will not be detected

					if (jmod < numRowsPixelMap)
					{
						pixelMap(jmod, column) += numExcessElectrons;

						// Make sure the pixel you move the excess electrons to
						// does not get saturated too

						if (pixelMap(jmod, column) > fullWellSaturationLimit)
						{
							numExcessElectrons = pixelMap(jmod, column) - fullWellSaturationLimit;
						}
					}
				}

				// Transfer excess electrons up

				jmod = row;
				numExcessElectrons = (pixelValue - fullWellSaturationLimit) / 2.0;    // ...and the rest of the excess electrons up

				while (numExcessElectrons > 0 && jmod >= 0)
				{
					pixelMap(jmod, column) -= numExcessElectrons;
					jmod--;

					// Electrons reaching the edge of the CCD will not be detected

					if (jmod >= 0)
					{
						pixelMap(jmod, column) += numExcessElectrons;

						// Make sure the pixel you move the excess electrons to does not get saturated too

						if (pixelMap(jmod, column) > fullWellSaturationLimit)
						{
							numExcessElectrons = pixelMap(jmod, column) - fullWellSaturationLimit;
						}
					}
				}
			}
		}
	}
}












/**
 * \brief: Apply the effect of the charge-transfer (in)efficiency to the
 *         pixel map. The serial register is assumed to have a CTE of 1, 
 *         unlike the CCD that has a CTE map.
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre Pixel unit in the smearing map: [electrons].
 * \pre No bias register map.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map.
 */
void Detector::applyCte()
{
	Log.debug("Detector: Applying charge transfer inefficiency");

	float cti = 1.0 - meanCte;

	// Computing the effects of CTE requires the use of a binomial distribution.
	// To speed up things, we first pre-compute some parts of this distribution.

	// Pre-compute the (natural) logarithms of the first N natural numbers

	vector<double> logs(numRowsPixelMap + subFieldZeroPointRow);
	iota(logs.begin(), logs.end(), 1.0);
	transform(logs.begin(), logs.end(), logs.begin(),
			ptr_fun<double, double>(log));

	// Compute the partial sums of these logarithms
	// sumOfLogsUpTo[i] contains log((i+1)!) = log(1) + ... + log(i+1)

	vector<double> sumOfLogsUpTo(numRowsPixelMap + subFieldZeroPointRow);
	partial_sum(logs.begin(), logs.end(), sumOfLogsUpTo.begin());

	arma::Row<float> readout;	// Readout strip

	// Loop over all rows in the pixel map (starting at the row farthest away from
	// the readout register)

	for (int row = numRowsPixelMap - 1; row >= 0; row--)
	{

		// Reset the readout register

		readout.zeros(numColumnsPixelMap);

		// Each row picks up flux that is left behind when transferring the rows
		// that are closer to the readout register, row-by-row to the readout
		// register (these rows are looped over via the "index" variable - note
		// that the detector zeropoint is added to it!).

		for (unsigned int index = subFieldZeroPointRow;
				index < row + subFieldZeroPointRow; index++)
		{
			const double factor1 = pow(meanCte, index + 1)
					* pow(cti, row + subFieldZeroPointRow - index);

			if ((row + subFieldZeroPointRow == 0) || (index == 0))
			{
				readout += pixelMap(index - subFieldZeroPointRow,
						arma::span::all) * factor1;
			} else
			{
				const double factor2 = exp(
						sumOfLogsUpTo[row + subFieldZeroPointRow - 1]
								- sumOfLogsUpTo[row + subFieldZeroPointRow
										- index - 1]
								- sumOfLogsUpTo[index - 1]);
				readout += pixelMap(index - subFieldZeroPointRow,
						arma::span::all) * factor1 * factor2;
			}
		}

		pixelMap(row, arma::span::all) = readout(0, arma::span::all);
	}

	// BELOW: OLD IMPLEMENTATION

//	// Create a map in which we will shift the rows of the pixel map one-by-one
//	// towards the readout register.  Bear in mind that the bottom row of the
//	// sub-field is not necessarily right next to the readout register (the
//	// distance between the two is subFieldZeroPointRow).
//
//	arma::Mat<float> shiftMap;
//	shiftMap.zeros(subFieldZeroPointRow + numRowsPixelMap, numColumnsPixelMap);
//	shiftMap.submat(arma::span(subFieldZeroPointRow, subFieldZeroPointRow + numRowsPixelMap - 1), arma::span::all) = pixelMap;
//
//	// The readout register
//
//	arma::Row<float> readoutStrip;
//	readoutStrip.zeros(numColumnsPixelMap);
//
//	// Array filled with ones (needed for the CTI)
//
//	arma::Row<float> ones;
//	ones.ones(numColumnsPixelMap);
//
//	// Shift all the rows down (i.e. towards the readout register) one-by-one
//	// Keep on doing this until all rows have been read out.
//
//	for (int shiftIndex = 0;
//			shiftIndex < numRowsPixelMap + subFieldZeroPointRow; shiftIndex++)
//	{
//
//		// Shift the bottom row to the readout strip
//
//		readoutStrip = cteMap(0, arma::span::all) * shiftMap(0, arma::span::all);
//
//		if (shiftIndex >= subFieldZeroPointRow)
//		{
//			pixelMap(shiftIndex - subFieldZeroPointRow, arma::span::all) = readoutStrip(0, arma::span::all);
//		}
//
//		// Shift all other rows one row down (i.e. closer to the readout register)
//
//		for (int row = 0; row < subFieldZeroPointRow + numRowsPixelMap - 1; row++)
//		{
//			shiftMap(row, arma::span::all) = (ones-cteMap(row, arma::span::all))
//					                           * shiftMap(row, arma::span::all)	// Left behind when shifting row down (CTI = 1 - CTE)
//					                         + cteMap(row + 1, arma::span::all)
//							                   * shiftMap(row + 1, arma::span::all);	// Transferred (CTE)
//		}
//	}
}











/**
 * \brief: Apply the effect of readout smearing to the pixel and the smearing map.
 *         This effect is due to the absence of a shutter (common in space-based 
 *         instruments) - the CCD still receives light during frame transfer.
 *         The flux of each pixel is affected by the flux of the pixels
 *         in the same column.  Because the CCD is exposed during the whole 
 *         readout and multiple exposures are created, also the pixels further 
 *         away from the readout register are affected.
 *
 * NOTES: A smearing map is created and will be used in photometry to remove 
 *        the smearing effect from the pixel map.
 *
 * \param exposureTime: Exposure time [s].
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre Pixel unit in the smearing map: [electrons].
 * \pre No bias register map.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map.
 */
void Detector::applyOpenShutterSmearing(float exposureTime)
{
	// Average out the fluxes in the pixel map per column and make sure it is
	// scaled with the readout time instead of with the exposure time.

	arma::Row<float> openShutterSmearing = arma::sum(pixelMap, 0);
	float factor = (readoutTime / exposureTime) / numRowsPixelMap;
	openShutterSmearing *= factor;

	// Add the effect of the open-shutter smearing to the pixel map

	for (unsigned int row = 0; row < numRowsPixelMap; row++)
	{
		pixelMap(row, arma::span::all) += openShutterSmearing;
	}

	// Add the effect of the open-shutter smearing to the smearing map

	for (unsigned int row = 0; row < numRowsSmearingMap; row++)
	{
		smearingMap(row, arma::span::all) += openShutterSmearing;
	}
}













/**
 * \brief: Apply the readout noise to the pixel map and initialises the
 *         bias map. 
 * 
 * NOTE: Readout noise occurs due to the imperfect nature of the CCD amplifiers.  
 *       When the electrons are transferred to the amplifier, the induced voltage
 *       is measured. However, this measurement is not perfect, but gives a value 
 *       which is on average correct, with the readout noise as standard deviation.
 *       So readout noise is a measure of this scatter around the true value.
 *       Its value is expressed in electrons as the packet of charge is made up of 
 *       electrons.
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre Pixel unit in the smearing map: [electrons].
 * \pre No bias register map.
 * \pre Readout noise expressed in [electrons].
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post Pixel unit in the bias register map: [electrons].
 * \post Initialised the bias register map with readout noise.
 */
void Detector::addReadoutNoise()
{
	Log.debug("Detector: adding readout noise");

	readoutNoiseDistribution = normal_distribution<double>(0.0, readoutNoise);

	// Add readout noise to the pixel map

	for (unsigned int row = 0; row < numRowsPixelMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsPixelMap; column++)
		{
			pixelMap(row, column) += readoutNoiseDistribution(readoutNoiseGenerator);
		}
	}

	// Initialise the bias register map with readout noise

	for (unsigned int row = 0; row < numRowsBiasMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsPixelMap; column++)
		{
			biasMap(row, column) += readoutNoiseDistribution(readoutNoiseGenerator);
		}
	}
}










/**
 * \brief: Divide the bias register, smearing, and pixel map by the detector gain.
 *         This converts these three maps from electrons to ADU:
 *
 * \pre Pixel unit in the pixel, smearing, and bias register maps: [electrons].
 *
 * \post Pixel unit in the pixel, smearing, and bias register maps: [ADU].
 */
void Detector::applyGain()
{
	Log.debug("Detector: applying gain");

	// Divide the pixel, bias register, and smearing map by the gain

	pixelMap /= gain;
	biasMap /= gain;
	smearingMap /= gain;
}











/**
 * \brief: Add the electronic offset (i.e. bias level) to the pixel map,
 *         smearing map, and bias map.
 *
 * \pre Pixel unit in the pixel, smearing, and bias maps: [ADU].
 * \pre Electronic offset (i.e. bias level) expressed in [ADU].
 *
 * \post Pixel unit in the pixel, smearing, and bias register maps: [ADU].
 * \post Electronic offset added to the pixel, smearing, and bias register maps.
 */

void Detector::addElectronicOffset()
{
	Log.debug("Detector: adding a bias");

	// Add the electronic offset to the pixel, bias register, and smearing maps

	pixelMap += electronicOffset;
	biasMap += electronicOffset;
	smearingMap += electronicOffset;
}












/**
 * \brief: Apply the effect of digital saturation to the pixel map,
 *         smearing map, and bias register map. This means that the pixel values in
 *         these maps (expressed in [ADU / pixel]) are topped off to the digital saturation
 *         limit of the detector (also expressed in [ADU / pixel]).
 *
 * \pre Pixel unit in the pixel, smearing, and bias maps: [ADU].
 * \pre Digital saturation limit expressed in [ADU / pixel].
 *
 * \post Pixel unit in the pixel, smearing, and bias register maps: [ADU].
 */
void Detector::applyDigitalSaturation()
{
	Log.debug("Detector: applying digital saturation");

	// Top off the values in the pixel map

	pixelMap(arma::find(pixelMap > digitalSaturationLimit)).fill(digitalSaturationLimit);

	// Top off the values in the bias register map

	biasMap(arma::find(biasMap > digitalSaturationLimit)).fill(digitalSaturationLimit);

	// Top off the values in the smearing map

	smearingMap(arma::find(smearingMap > digitalSaturationLimit)).fill(digitalSaturationLimit);
}












/**
 * \brief Compute the planar (x,y) coordinates in the FP' reference system (not the FP system) 
 *        given the (real-valued) pixel coordinates on the CCD.
 *        
 * \param row     Row coordinate, real-valued (e.g. 3.5)    [pix]
 * \param column  Column coordinate, real-valued (e.g. 8.3) [pix]
 * 
 * \return (xFPprime, yFPprime)  A pair of (x,y) coordinates in the FP' reference system [mm]
 */

pair<double, double> Detector::pixelToPlanarFocalPlaneCoordinates(double row, double column)
{
    // Convert the pixel coordinates into [mm] coordinates
    // The pixelSize is expressed in [micron].

    double xCCDmm = row * pixelSize / 1000.0;
    double yCCDmm = column * pixelSize / 1000.0;

    // Convert the CCD coordinates into FP' coordinates [mm]
    // Note: orientationAngle is in [rad], originOffsetX and originOffsetY in mm

    double xFPprime = originOffsetX + xCCDmm * cos(orientationAngle) - yCCDmm * sin(orientationAngle);
    double yFPprime = originOffsetY + xCCDmm * sin(orientationAngle) + yCCDmm * cos(orientationAngle);

    // That's it

    return make_pair(xFPprime, yFPprime);
}










/**
 * \brief Compute the (real-valued) pixel coordinates of the star on the CCD, given the 
 *        planar (x,y) coordinates in the FP' reference system (not the FP system)
 *
 * \param xFPprime  planar x-coordinate of the point in the FP' reference system  [mm]
 * \param yFPprime  planar y-coordinate of the point in the FP' reference system  [mm]
 * 
 * \return (row, column)  Row and column pixel coordinates of the point (real-valued) [pix]
 */

pair<double, double> Detector::planarFocalPlaneToPixelCoordinates(double xFPprime, double yFPprime)
{
	// Convert the FP' coordinates into CCD coordinates [mm]

    double xCCDmm =  (xFPprime-originOffsetX) * cos(orientationAngle) + (yFPprime-originOffsetY) * sin(orientationAngle);
    double yCCDmm = -(xFPprime-originOffsetX) * sin(orientationAngle) + (yFPprime-originOffsetY) * cos(orientationAngle);

    // Convert the [mm] coordinates into pixel coordinates
    // Note: the pixel size is expressed in [micron]

    double column = xCCDmm / pixelSize * 1000.0;
    double row = yCCDmm / pixelSize * 1000.0;

    // That's it

    return make_pair(row, column);
}












/**
 * \brief  Return the focal plane coordinates of the center pixel of the subfield
 * 
 * \return (x,y)   focal plane coordinates in the FP' reference system [mm]
 */

pair<double, double> Detector::getPlanarFocalPlaneCoordinatesOfSubfieldCenter()
{
	double centerRow = subFieldZeroPointRow + numRowsPixelMap / 2.0;
	double centerCol = subFieldZeroPointColumn + numColumnsPixelMap / 2.0;

	return pixelToPlanarFocalPlaneCoordinates(centerRow, centerCol);
}












/**
 * \brief  Return the distance between the two diagonally opposite corner points of the
 *         subfield in [mm] on the focal plane.
 * 
 * \return length  Diagonal distance between the two opposite corners [mm]
 */

double Detector::getDiagonalLengthOfSubfield()
{
	// Define the pixel row and column coordinates of two corners diagonally opposite to each other

	double corner1Row = subFieldZeroPointRow + numRowsPixelMap;
	double corner1Col = subFieldZeroPointColumn;

	double corner2Row = subFieldZeroPointRow;
	double corner2Col = subFieldZeroPointColumn + numColumnsPixelMap;

	// Compute their corresponding (x,y) coordinates in the focal plane FP'

	double corner1X, corner1Y;
	double corner2X, corner2Y;

	tie(corner1X, corner1Y) = pixelToPlanarFocalPlaneCoordinates(corner1Row, corner1Col);
	tie(corner2X, corner2Y) = pixelToPlanarFocalPlaneCoordinates(corner2Row, corner2Col);

	// The diagonal length is simply the distance between the two corner points
	// computed using Pythagoras

	double diagonalLength = sqrt(pow(corner1X - corner2X, 2) + pow(corner1Y - corner2Y, 2));

	return diagonalLength;
}












/**
 * \brief: Creates the group(s) in the HDF5 file where the detector specific
 *         information will be stored.  These groups have to be created once,
 *         at the very beginning.
 */
void Detector::initHDF5Groups()
{
	Log.debug("Detector: initialising HDF5 groups");

	hdf5File.createGroup("/Images");
	hdf5File.createGroup("/BiasMaps");
	hdf5File.createGroup("/SmearingMaps");
	hdf5File.createGroup("/Flatfield");
}












/**
 * \brief: Writes the pixel map for the HDF5 file.
 */
void Detector::writePixelMapToHDF5()
{
	stringstream myStream;
    myStream << "image" << setfill('0') << setw(6) << imageNr;
    string imageName = myStream.str();

    // Add the image to the "Images" group

    hdf5File.writeArray("/Images", imageName, pixelMap);

    // Increment the counter for the next image

    imageNr++;
}
