/*
 * subfield.cpp
 *
 *  Created on: 14 Dec 2015
 *      Author: sara
 */

# include "subfield.h"

/**
 * Constructor.
 *
 * Stores all information to convert from focal plane coordinates [mm] to
 * (sub-)pixel coordinates in the sub-pixel map.
 *
 * @param configurationParameters: Configuration parameters for the sub-field.
 * @type configurationParameters: ConfigurationParameters
 *
 * @param detectorOriginOffsetX: Detector origin offset in the x-direction [mm].
 * @type detectorOriginOffsetX: double
 *
 * @param detectorOriginOffsetY: Detector origin offset in the y-direction [mm].
 * @type detectorOriginOffsetY: double
 *
 * @param detectorOrientation: Orientation angle of the detector w.r.t. the
 *        focal plane [degree].
 * @type detectorOrientation: double
 *
 * @param pixelSize: Detector pixel size [mm / pixel].
 * @type pixelSize: double
 */
SubField::SubField(ConfigurationParameters configurationParameters,
		double detectorOriginOffsetX, double detectorOriginOffsetY,
		double detectorOrientation, double pixelSize)
{

	// Number of sub-pixels per pixel

	// Detector origin offset in both directions [mm]

	this->detectorOriginOffsetX = detectorOriginOffsetX;
	this->detectorOriginOffsetY = detectorOriginOffsetY;

	// Pixel size [mm / pixel]

	this->pixelSize = pixelSize;

	// Detector orientation w.r.t. the focal plane [radians]

	this->detectorOrientation = detectorOrientation * M_PI / 180.0;

	// Number of pixels to add to each side of the sub-field to account for the
	// edge effect

//	this->numEdgePixels=...

	// Dimensions of the sub-pixel map, incl. the edge pixels

	this->subPixelMapSizeX = (this->subFieldSizeX + 2 * this->numEdgePixels)
			* this->numSubPixelsPerPixel;
	this->subPixelMapSizeY = (this->subFieldSizeY + 2 * this->numEdgePixels)
			* this->numSubPixelsPerPixel;

	this->subPixelMap = new double*[this->subPixelMapSizeY];

	for (unsigned int row = 0; row < this->subPixelMapSizeY; row++)
	{
		this->subPixelMap[row] = new double[this->subPixelMapSizeX];
	}
}

/**
 * Destructor.
 */
SubField::~SubField()
{

	// Destroy the sub-pixel map

	for (unsigned int row = 0; row < this->subPixelMapSizeY; row++)
	{
		delete[] this->subPixelMap[row];
	}

	delete[] this->subPixelMap;
}

/**
 * Method that resets the sub-pixel map, including the edge pixels that are
 * added on each side to account for the edge effect.
 *
 * @pre The sub-pixel map either does not exist or contains values that must be
 *      set to zero.
 *
 * @post The sub-pixel maps is either initialised or its values are set to zero.
 */
void SubField::reset()
{

	// Initialise sub-pixel map
	// Dimensions: this-> subPixelMapSize
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
void SubField::addFlux(double xCoords, double yCoords, double flux)
{

	// Detector origin offset (pixel level)

	double xOffset = (xCoords - this->detectorOriginOffsetX) / this->pixelSize;
	double yOffset = (yCoords - this->detectorOriginOffsetY) / this->pixelSize;

	// Detector orientation (pixel level)

	double column = xOffset * cos(this->detectorOrientation)
			- yOffset * sin(this->detectorOrientation);
	double row = xOffset * sin(this->detectorOrientation)
			+ yOffset * cos(this->detectorOrientation);

	// Sub-field incl. edge pixels (also correct for sub-field zeropoint)

	column = (column - this->subFieldZeroPointX + this->numEdgePixels)
			* this->numSubPixelsPerPixel;
	row = (row - this->subFieldZeroPointY + this->numEdgePixels)
			* this->numSubPixelsPerPixel;

	// Add flux in this->subPixelMap at (row, column)

	if (this->isInSubPixelMap(row, column))
	{
		this->subPixelMap[(int) row][(int) column] += flux;
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
bool SubField::isInSubPixelMap(double row, double column)
{
	return (column >= 0) && (row >= 0) && (column < this->subPixelMapSizeX)
			&& (row < this->subPixelMapSizeY);
}

/**
 * Method that adds the given flux value to (all sub-pixels of) the sub-pixel map.
 *
 * @param flux: Flux to add to the sub-pixel map.
 * @type flux: double
 */
void SubField::addFlux(double flux)
{
	for (unsigned int row = 0; row < this->subPixelMapSizeY; row++)
	{
		for (unsigned int column = 0; column < this->subPixelMapSizeX; column++)
		{

			this->subPixelMap[row][column] += flux;

		}
	}
}

/**
 * Method that convolves the sub-field with the given PSF.
 *
 * @param psf: PSF to convolve the sub-field map with.
 * @type psf: double**
 */
void SubField::convolveWithPsf(double **psf)
{

}

/**
 * Method that multiplies the sub-pixel map with the given array.  If the dimensions
 * do not match (the sub-pixel map contains edge pixels and the given array doesn't),
 * the central part of the sub-pixel map is multiplied with the given array.
 *
 * @param array: Array by with to multiply the sub-pixel map.
 * @type array: double**
 */
void SubField::multiply(double **array)
{

}

/**
 * Method that rebins the sub-pixel map to pixel level and crops the edge pixels
 * that were added on each side to account for the edge effect.
 */
void SubField::rebin()
{

}

