#include "detector.h"

/**
 * Constructor.  Creates a detector object, based on the given configuration
 * parameters and attaches it to the given camera.
 *
 * @param configurationParameters: Configuration parameters for the detector.
 * @type configurationParameters: ConfigurationParameters
 *
 * @param camera: Camera to which to attach the detector.
 * @type camera: Camera
 */
Detector::Detector(ConfigurationParameters configurationParameters,
		Camera camera)
{

	// Parse the parameters from the configuration file

	// Associate the camera

	// Initialise flatfield map

	// this->setFlatfieldMap(peak2PeakNoise, subPixelNoise, intraPixelWidth);

	// Initialise CTE map
}





/**
 * Destructor.
 *
 * @pre Sub-pixel map at pixel level, excl. edge pixels.
 * @pre CTE map at pixel level, excl. edge pixels.
 * @pre Flatfield map at sub-pixel level, excl. edge pixels.
 * @pre Bias register map at pixel level, excl. edge pixels.
 * @pre Smearing map at pixel level, excl. edge pixels.
 *
 * @post De-allocated memory of sub-field.
 * @post De-allocated memory of CTE map.
 * @post De-allocated memory of flatfield map.
 * @post De-allocated memory of bias register map.
 * @post De-allocated memory of smearing map.
 */
Detector::~Detector()
{

	// Destroy the sub-field and the CTE map

	for(unsigned int row = 0; row < subFieldSizeY; row++)
	{
		delete[] subField[row];
		delete cteMap[row];
	}

	delete[] subField;
	delete[] cteMap;

	// Destroy the flatfield map

	for(unsigned int row = 0; row < subFieldSizeY*numSubPixelsPerPixel; row++)
	{
		delete[] flatfieldMap[row];
	}

	delete[] flatfieldMap;

	// Delete the bias register map

	for(unsigned int row = 0; row < numBiasPrescanRows; row++)
	{
		delete[] biasRegisterMap[row];
	}

	delete[] biasRegisterMap;

	// Delete the smearing map

	for(unsigned int row = 0; row < numSmearingOverscanRows; row++)
	{
		delete[] smearingMap[row];
	}

	delete[] smearingMap;
}





/**
 * Method that resets the sub-field, the bias register map, and the smearing map.
 *
 * @pre Sub-field has either not been initialised yet, or is defined at pixel
 *      level and possibly contains values from previous exposure.
 * @pre Bias register map has either not been initialised yet, or has the correct
 *      dimensions but is possibly filled with values from previous exposure.
 * @pre Smearing map has either not been initialised yet, or has the correct
 *      dimensions but is possibly filled with values from previous exposure.
 *
 * @post Sub-field of correct dimensions (i.e. at sub-pixel level, incl. edge
 *       pixels), filled with zeroes.
 * @post Bias register map of correct dimensions, filled with zeroes.
 * @post Smearing map of correct dimensions, filled with zeroes.
 */
void Detector::reset()
{

	// Reset sub-field

	this->resetSubField();

	// Reset bias register map

	this->resetBiasRegisterMap();

	// Reset smearing map

	this->resetSmearingMap();

}





/**
 * Method that resets the sub-field.
 *
 * @pre Sub-field has either not been initialised yet, or is defined at pixel
 *      level and possibly contains values from previous exposure.
 *
 * @post Sub-field of correct dimensions (i.e. at sub-pixel level, incl. edge
 *       pixels), filled with zeroes.
 */
void Detector::resetSubField()
{
	// De-allocate memory (pixel level, excl. edge pixels)

	if (subField != nullptr)
	{
		for (unsigned int row = 0; row < subFieldSizeY; row++)
		{
			delete[] subField[row];
		}

		delete[] subField;
	}

	// Allocate memory (sub-pixel level, incl. edge pixels)

	subField = new double*[subFieldSizeY];

	for (unsigned int row = 0; row < subPixelMapSizeY; row++)
	{
		subField[row] = new double[subFieldSizeX];
	}


	// Set all values to zero

	memset(subField, 0,
			sizeof(subField[0][0]) * subPixelMapSizeX * subPixelMapSizeY);
}





/**
 * Mehod that resets the bias register map, i.e. all values are set to zero.
 *
 * @pre Bias register map has either not been initialised yet, or has the correct
 *      dimensions but is possibly filled with values from previous exposure.
 *
 * @post Bias register map of correct dimensions, filled with zeroes.
 */
void Detector::resetBiasRegisterMap()
{
	// Allocate memory

	if(biasRegisterMap == nullptr)
	{
		biasRegisterMap = new double*[numBiasPrescanRows];

		for(unsigned int row = 0; row < numBiasPrescanRows; row++)
		{
			biasRegisterMap[row] = new double[subFieldSizeX];
		}
	}

	// Set all values to zero

	memset(biasRegisterMap, 0,
			sizeof(biasRegisterMap[0][0]) * numBiasPrescanRows * subFieldSizeX);
}





/**
 * Method that resets the smearing map, i.e. all values are set to zero.
 *
 * @pre Smearing map has either not been initialised yet, or has the correct
 *      dimensions but is possibly filled with values from previous exposure.
 *
 * @post Smearing map of correct dimensions, filled with zeroes.
 */
void Detector::resetSmearingMap()
{
	// Allocate memory

	if (smearingMap == nullptr)
	{
		smearingMap = new double*[numSmearingOverscanRows];

		for (unsigned int row = 0; row < numSmearingOverscanRows; row++)
		{
			smearingMap[row] = new double[subFieldSizeX];
		}
	}

	// Set all values to zero

	memset(smearingMap, 0,
			sizeof(smearingMap[0][0]) * numSmearingOverscanRows
					* subFieldSizeX);
}





/**
 * Method that takes an exposure with the detector starting at the given time.
 *
 * The light is integrated during the given exposure time, during which the
 * detector suffers from the effects of jitter and telescope drift (i.e. thermo-
 * elastic variations).  The background is assumed uniform for the whole
 * sub-field.
 *
 * Afterwards, the collected light is read out, convolving the image with the
 * PSF of the camera and adding various noise effects.
 *
 * @param startTime: Starting time of the exposure [s].
 * @type startTime: double
 *
 * @param exposureTime: Duration of the exposure [s].
 * @type exposureTime: double
 */
void Detector::takeExposure(double startTime, double exposureTime)
{

	// Integration (jitter + drift) + background contribution

	this->integrateLight(startTime, exposureTime);

	// Readout: convolution with PSF + noise effects

	this->readOut();
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
 * @param startTime: Starting time of the exposure for which jitter must be
 *        applied [s].
 *
 * @pre No sub-pixel map
 * @pre No pixel map
 * @pre No smearing map
 * @pre No bias register map
 *
 * @post Pixel unit of the sub-pixel map: [e-]
 * @post Re-calculated star position for each jitter step and summed the flux
 *       in the sub-pixel map, which introduces blurring. Also the contribution
 *       of the background (zodiacal + galactic) is added to the sub-pixel map.
 * @post No pixel map
 * @post No smearing map
 * @post No bias register map
 */
void Detector::integrateLight(double startTime, double exposureTime)
{

	// Reset the sub-field (i.e. get rid of the previous exposure, by zeroing the entire sub-field)

	this->reset();

	// Integration (incl. jitter) + background

	this->camera.exposeSubField(this);

	// Apply flatfield (at sub-pixel level)

	this->applyFlatfield();

	// Rebin

	this->rebin();
}





/**
 * Method that adds the given flux value to the value of the sub-pixel that
 * corresponds to the given coordinates in the focal plane.
 *
 * @param xCoords: The x-coordinate of the sub-pixel in the focal plane [mm].
 * @type xCoords: double
 *
 * @param yCoords: The y-coordinate of the sub-pixel in the focal plane [mm].
 * @type yCoords: double
 *
 * @param flux: Flux to add to the sub-pixel map.
 * @type flux: double
 */
void Detector::addFlux(double xCoords, double yCoords, double flux)
{

	// Detector origin offset (pixel level)

	double xOffset = (xCoords - originOffsetX) / pixelSize;
	double yOffset = (yCoords - originOffsetY) / pixelSize;

	// Detector orientation (pixel level)

	double column = xOffset * cos(orientationAngle)
			- yOffset * sin(orientationAngle);

	double row = xOffset * sin(orientationAngle)
			+ yOffset * cos(orientationAngle);

	// Sub-field incl. edge pixels (also correct for sub-field zeropoint)

	column = (column - subFieldZeroPointX + numEdgePixels)
			* numSubPixelsPerPixel;

	row = (row - subFieldZeroPointY + numEdgePixels) * numSubPixelsPerPixel;

	// Add flux in this->subPixelMap at (row, column)

	if (this->isInSubField(row, column))
	{
		subField[(int) row][(int) column] += flux;
	}
}





/**
 * Method that checks whether the given (row, column) coordinates are in the
 * sub-pixel map.
 *
 * @param row: Row coordinate.
 * @type row: double
 *
 * @param column: Column coordinate.
 * @type column: double
 *
 * @return True if the given (row, column) coordinates are in the sub-pixel map;
 *         false otherwise.
 */
bool Detector::isInSubField(double row, double column)
{

	return (column >= 0) && (row >= 0) && (column < subPixelMapSizeX)
			&& (row < subPixelMapSizeY);

}





/**
 * Method that adds the given flux value to (all sub-pixels of) the sub-pixel map.
 *
 * @param flux: Flux to add to the sub-pixel map.
 * @type flux: double
 */
void Detector::addFlux(double flux)
{
	for (unsigned int row = 0; row < subPixelMapSizeY; row++)
	{
		for (unsigned int column = 0; column < subPixelMapSizeX; column++)
		{

			subField[row][column] += flux;

		}
	}
}





void Detector::applyFlatfield()
{

}





/**
 * Method that rebins the sub-pixel map to pixel level and crops the edge pixels
 * that were added on each side to account for the edge effect.
 */
void Detector::rebin()
{
	// Rebinning is only required of the pixels are divided into sub-pixels

	if (numSubPixelsPerPixel > 1)
	{

		// Allocate memory

		double **pixelMap = new double*[subFieldSizeY];

		for (unsigned row = 0; row < subFieldSizeY; row++)
		{
			pixelMap[row] = new double[subFieldSizeX];
		}

		for (unsigned int row = 0; row < subFieldSizeY; row++)
		{
			for (unsigned int column = 0; column < subFieldSizeX; column)
			{

				double sum = 0;

				for (unsigned int r = row * numSubPixelsPerPixel;
						r < (row + 1) * numSubPixelsPerPixel; r++)
				{
					for (unsigned int c = column * numSubPixelsPerPixel;
							c < (column + 1) * numSubPixelsPerPixel; c++)
					{
						sum += subField[r][c];
					}
				}

				pixelMap[row][column] = sum;
			}
		}

		for(unsigned int row = 0; row < subPixelMapSizeY; row++)
		{
			delete[] subField[row];
		}

		delete[] subField;

		subField = pixelMap;
		pixelMap = nullptr;
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
 */
void Detector::readOut()
{

	// Apply quantum efficiency
	// Pixel units before: [photons]
	// Pixel units after: [electrons]

	this->applyQuantumEfficiency();

	// Apply poisson distributed photon noise
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	this->addPhotonNoise();

	// Apply full well saturation. A pixel has a maximum capacity of electrons (the full well capacity).
	// If photons free more electrons, the pixel saturates, and the electrons flow in the pixels above and below in
	// the same column (potential barriers are smallest in that direction).
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	this->applyFullWellSaturation();

	// Simulate the effects of the Charge Transfer Inefficiency (CTI). When the
	// CCD is read out, row after row, a part of the charge is always left behind
	// which then dribbles into the trailing pixels. This causes each star to have
	// a small "tail". Only visible when the CTI = 1 - CTE is poor.
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	this->applyCte();

	// Apply the effects of readout smearing due to an open shutter. Because there is no shutter,
	// the pixels are still receiving photons from the sky, while they are being transfered towards
	// the readout register.

	this->applyOpenShutterSmearing();

	// Each time the amplifier reads out a pixel, a tiny bit of noise is added.
	// Add the readout noise.
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	this->addReadoutNoise();

	// Apply the gain, to increase the dynamic range of the detector.
	// Pixel units before: [electrons]
	// Pixel units after: [ADU]

	this->applyGain();

	// Take into account the bias level. I.e. add the constant "zero" level
	// introduced by the amplifier.
	// Pixel units before: [ADU]
	// Pixel units after: [ADU]

	this->addElectronicOffset();

	// Take into acount digital saturation. If even after dividing by the gain
	// the number of ADUs in a pixel is still higher than the analog-digital
	// converter (ADC) can represent with its fixed amount of bits, clip all
	// values that are too high to the saturation level of the ADC.
	// Pixel units before: [ADU]
	// Pixel units after: [ADU]

	this->applyDigitalSaturation();
}





void Detector::applyQuantumEfficiency()
{

}

/**
 * Method that adds photon noise (i.e. shot noise) to the pixel map.  This type
 * of noise occurs because of the discrete/quantised nature of the electric
 * charge carried by the electrons (when counting them as representatives of
 * photons hitting the detectors).  It follows a Poisson distribution and each
 * pixel is treated independent of the other pixels.
 *
 * @pre Pixel unit in the pixel map: [e-]
 * @pre No smearing map
 * @pre No bias register map
 *
 * @post Pixel unit in the pixel map: [e-]
 * @post Added photon noise to the pixel map
 * @post Pixel unit in the smearing map: [e-]
 * @post Added photon noise to the smearing map
 * @post No bias register map
 */
void Detector::addPhotonNoise()
{

	// Add photon noise to the pixel map

	// Iniitialise the smearing map with photon noise
}





/**
 * Method that applies the effect of full-well saturation (i.e. blooming) to the
 * pixel map.  If a pixel receives more electrons than the full-well saturation
 * limit (expressed in [e-/pixel]), the additional electrons flow evenly
 * distributed in positive and negative charge-transfer direction.  Electrons
 * reaching the edge of the CCD will not be detected.
 *
 * @pre Pixel unit in the pixel map: [e-]
 * @pre Pixel unit in the smearing map: [e-]
 * @pre No bias register map
 * @pre Full-well saturation limit expressed in [e-]
 *
 * @post Pixel unit in the pixel map: [e-]
 * @post Effect of full-well saturation (i.e. blooming) applied to the pixel map
 * @post Pixel unit in the smearing map: [e-]
 * @post No bias register map
 */
void Detector::applyFullWellSaturation()
{

}





/**
 * Method that applies the effect of the charge-transfer (in)efficiency to the
 * pixel map. Assumed is that the serial register has a CTE of 1, unlike the
 * CCD that have a CTE map.
 *
 * @pre Pixel unit in the pixel map: [e-]
 * @pre Pixel unit in the smearing map: [e-]
 * @pre No bias register map
 *
 * @post Pixel unit in the pixel map: [e-]
 * @post Applied the effect of the charge-transfer (in)efficiency of the CCD to
 *       the pixel map.
 * @post Pixel unit in the smearing map: [e-]
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
 * @pre Pixel unit in the pixel map: [e-]
 * @pre Pixel unit in the smearing map: [e-]
 * @pre No bias register map
 *
 * @post Pixel unit in the pixel map: [e-]
 * @post Applied the effect of readout smearing to the pixel map
 * @post Pixel unit in the smearing map: [e-]
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
 * Its value is expressed in electrons, because, after all, the packet of charge
 * is made up of electrons.
 *
 * @pre Pixel unit in the pixel map: [e-]
 * @pre Pixel unit in the smearing map: [e-]
 * @pre No bias register map
 * @pre Readout noise expressed in [e-]
 *
 * @post Pixel unit in the pixel map: [e-]
 * @post Added readout noise to the pixel map
 * @post Pixel unit in the smearing map: [e-]
 * @post Pixel unit in the bias register map: [e-]
 * @post Initialised the bias register map with readout noise
 */
void Detector::addReadoutNoise()
{

	// Normal<double, ranlib::MersenneTwister, ranlib::independentState> randomMap(
	//		0, this->getReadoutNoise());
	// randomMap.seed(p_DataSet->getSeedReadOutNoise());

	// add randomMap.random() to all pixel values!!

	// Add readout noise to the pixel map

	// Initialise the bias register map with readout noise
}





/**
 * Method that divides the bias register map, smearing map, and pixel map by the
 * detector gain. This relates the number of electrons per pixels to the number
 * of counts (i.e. ADU) per pixel and converts these three maps from electrons
 * to ADU:
 * <ul>
 * <li>bias register map [ADU / pixel] = bias register map [e- / pixel] / gain [e- / ADU]</li>
 * <li>smearing map [ADU / pixel] = smearing map [e- / pixel] / gain [e- / ADU]</li>
 * <li>pixel map [ADU / pixel] = pixel map [e- / pixel] / gain [e- / ADU]</li>
 * </ul>
 *
 * @pre Pixel unit in the pixel map: [e-]
 * @pre Pixel unit in the smearing map: [e-]
 * @pre Pixel unit in the bias register map: [e-]
 *
 * @post Pixel unit in the pixel map: [ADU]
 * @post Pixel unit in the smearing map: [ADU]
 * @post Pixel unit in the bias register map: [ADU]
 */
void Detector::applyGain()
{

	// Multiply the pixel map with the gain

	// Multiply the smearing map with the gain

	// Multiply the bias register map with the gain
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

	// Add the electronic offset to the smearing map

	// Add the electronic offset to the bias register map

}





/**
 * Method that applies the effect of digital saturation to the pixel map,
 * smearing map, and bias register map.  This means that the pixel values in
 * these maps (expressed in [ADU/pixel]) are topped off to the digital saturation
 * limit of the detector (also expressed in [ADU/pixel]).
 *
 * @pre Pixel unit in the pixel map: [ADU]
 * @pre Pixel unit in the smearing map: [ADU]
 * @pre Pixel unit in the bias register map: [ADU]
 * @pre Digital saturation limit expressed in [ADU/pixel]
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

	// Top off the values in the smearing map

	// Top off the values in the bias register map

}





///**
// * Method that returns the sub-field of the detector that is modelled in more
// * detail.
// *
// * @return Sub-field of the detector that is modelled in more detail.
// * @rtype SubField
// */
//SubField Detector::getSubField()
//{
//
//}
//
//
//
//
//
///**
// * Method that returns the offset of the detector origin from the centre of the
// * optical plane (i.e. optical axis) [mm].
// *
// * @return Offset of the detector origin from the centre of the optical plane [mm].
// * @rtype double
// */
//double Detector::getOriginOffsetX()
//{
//
//	return this->originOffsetX;
//
//}
//
//
//
//
//
///**
// * Method that returns the offset of the detector origin from the centre of the
// * optical plane (i.e. optical axis) [mm].
// *
// * @return Offset of the detector origin from the centre of the optical plane [mm].
// * @rtype double
// */
//double Detector::getOriginOffsetY()
//{
//
//	return this->originOffsetY;
//
//}
//
//
//
//
//
///**
// * Method that returns the pixel size of the detector [mm / pixel].
// *
// * @return Pixel size of the detector [mm / pixel].
// * @rtype unsigned double
// */
//double Detector::getPixelSize()
//{
//
//	return this->pixelSize;
//
//}
//
//
//
//
//
///**
// * Method that returns the orientation angle of the detector on the sky, measured
// * counterclockwise [degrees].
// *
// * @return Orientation angle of the detector on the sky [degrees].
// * @rtype double
// */
//double Detector::getOrientationAngle()
//{
//
//	return this->orientationAngle;
//
//}
