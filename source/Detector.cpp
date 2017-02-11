#include "Detector.h"

/**
 * \brief Constructor.
 * 
 * \details
 * 
 * The constructor initializes the groups in the HDF5 file where the different maps (i.e. pixel map,
 * bias register map, smearing map, etc.) will be saved. 
 * 
 * The following maps are initialized to zero:
 * 
 * pixelMap 
 * subPixelMap
 * biasMap
 * smearingMap
 * flatfieldMap
 * throughputMap
 * cteMap
 * 
 * The flatfieldMap is filled at sub-pixel level, the throughputMap and cteMap are filled at pixel level.
 * 
 * \param configParam    Configuration parameters for the detector.
 * \param hdf5file       HFD5 file to write the detector images to.
 * \param camera         Camera to which to attach the detector.
 */

Detector::Detector(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera)
: HDF5Writer(hdf5file), 
  includeFlatfield(true), 
  includePhotonNoise(true), 
  includeReadoutNoise(true),
  includeCTIeffects(true), 
  includeOpenShutterSmearing(true), 
  includeQuantumEfficiency(true),
  includeVignetting(true),
  includeParticulateContamination(true),
  includeMolecularContamination(true),
  writeSubPixelImagesToHDF5(false),
  includeFullWellSaturation(true),
  includeDigitalSaturation(true),
  psfWasSet(false), 
  internalTime(0.0), camera(camera), imageNr(0), subPixelImageNr(0)
{
	// Parse the parameters from the configuration file.

	configure(configParam);

	// Create the groups in the HDF5 file where the different maps (i.e. pixel map,
	// bias register map, smearing map, etc.) will be saved. This needs to be done
	// BEFORE other methods write arrays to HDF5.

	initHDF5Groups();

	// Allocate memory for the different maps

	pixelMap.zeros(numRowsPixelMap, numColumnsPixelMap);
	subPixelMap.zeros(numRowsSubPixelMap, numColumnsSubPixelMap);
	biasMap.zeros(numRowsBiasMap, numColumnsPixelMap);
	smearingMap.zeros(numRowsSmearingMap, numColumnsPixelMap);
	flatfieldMap.ones(numRowsSubPixelMap, numColumnsSubPixelMap);
	throughputMap.ones(numRowsPixelMap, numColumnsPixelMap);

	// Generate the flatfield map 

	generateFlatfieldMap();

	// Generate the throughput map (vignetting + polarisation + particulate & molecular contamination + quantum efficiency)

	generateThroughputMap();


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
	flushOutput();
}











/**
 * \brief Configure the Detector object using the ConfigurationParameters
 * 
 * \param configParam: the configuration parameters 
 **/

 void Detector::configure(ConfigurationParameters &configParam)
 {
 	// Configuration parameters for the CCD detector

    originOffsetX              			= configParam.getDouble("CCD/OriginOffsetX");
    originOffsetY              			= configParam.getDouble("CCD/OriginOffsetY");
    orientationAngle           			= deg2rad(configParam.getDouble("CCD/Orientation"));
    numRows                    			= configParam.getInteger("CCD/NumRows");
    numColumns                 			= configParam.getInteger("CCD/NumColumns");
    pixelSize                  			= configParam.getDouble("CCD/PixelSize");
    gain                       			= configParam.getInteger("CCD/Gain");
    quantumEfficiency          			= configParam.getDouble("CCD/QuantumEfficiency/Efficiency");
    refAngleQuantumEfficiency            = configParam.getDouble("CCD/QuantumEfficiency/RefAngle");
    expectedValueQuantumEfficiency       = configParam.getDouble("CCD/QuantumEfficiency/ExpectedValue");
    fullWellSaturationLimit    			= configParam.getLong("CCD/FullWellSaturation");
    digitalSaturationLimit     			= configParam.getLong("CCD/DigitalSaturation");
    readoutNoise              			= configParam.getDouble("CCD/ReadoutNoise");
    electronicOffset           			= configParam.getInteger("CCD/ElectronicOffset");
    readoutTime                			= configParam.getDouble("CCD/ReadoutTime");
    flatfieldNoiseAmplitude    			= configParam.getDouble("CCD/FlatfieldPtPNoise");
    expectedValueVignetting              = configParam.getDouble("CCD/Vignetting/ExpectedValue");
    particulateContaminationEfficiency 	= configParam.getDouble("CCD/Contamination/ParticulateContaminationEfficiency");
    molecularContaminationEfficiency     = configParam.getDouble("CCD/Contamination/MolecularContaminationEfficiency");

    CTImodel                   = configParam.getString("CCD/CTI/Model");
    if (CTImodel == "Simple")
    {
        meanCte                = configParam.getDouble("CCD/CTI/Simple/MeanCTE");
    }
    else if (CTImodel == "Short2013")
    {
        beta                    = configParam.getDouble("CCD/CTI/Short2013/Beta");
        temperature             = configParam.getDouble("CCD/CTI/Short2013/Temperature");
        numTrapSpecies          = configParam.getInteger("CCD/CTI/Short2013/NumTrapSpecies");   
        trapDensity             = configParam.getDoubleVector("CCD/CTI/Short2013/TrapDensity");
        trapCaptureCrossSection = configParam.getDoubleVector("CCD/CTI/Short2013/TrapCaptureCrossSection");
        releaseTime             = configParam.getDoubleVector("CCD/CTI/Short2013/ReleaseTime");
    }
    else
    {
        Log.error("Detector::configure(): Unkown CTI model specification in configuration file: "  + CTImodel);
        throw ConfigurationException("Detector: Unkown CTI model specification in configuration file");
    }

    polarizationEfficiency = configParam.getDouble("CCD/Polarization/Efficiency");
    refAnglePolarization = configParam.getDouble("CCD/Polarization/RefAngle");
    expectedValuePolarization = configParam.getDouble("CCD/Polarization/ExpectedValue");



    includeFlatfield                = configParam.getBoolean("CCD/IncludeFlatfield");
    includeParticulateContamination = configParam.getBoolean("CCD/IncludeParticulateContamination");
    includeMolecularContamination   = configParam.getBoolean("CCD/IncludeMolecularContamination");
    includePhotonNoise              = configParam.getBoolean("CCD/IncludePhotonNoise");
    includeReadoutNoise             = configParam.getBoolean("CCD/IncludeReadoutNoise");
    includeCTIeffects               = configParam.getBoolean("CCD/IncludeCTIeffects");
    includeOpenShutterSmearing      = configParam.getBoolean("CCD/IncludeOpenShutterSmearing");
    includeQuantumEfficiency        = configParam.getBoolean("CCD/IncludeQuantumEfficiency");
    includeVignetting               = configParam.getBoolean("CCD/IncludeVignetting");
    includePolarization             = configParam.getBoolean("CCD/IncludePolarization");
    writeSubPixelImagesToHDF5       = configParam.getBoolean("CCD/WriteSubPixelImagesToHDF5");
    includeConvolution              = configParam.getBoolean("CCD/IncludeConvolution");
    includeFullWellSaturation       = configParam.getBoolean("CCD/IncludeFullWellSaturation");
    includeDigitalSaturation        = configParam.getBoolean("CCD/IncludeDigitalSaturation");
    writeSubPixelImagesToHDF5       = configParam.getBoolean("CCD/WriteSubPixelImagesToHDF5");

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

	// Derive the dimensions of the sub-pixel map

	numRowsSubPixelMap    = numRowsPixelMap    * numSubPixelsPerPixel;	// TODO Add edge pixels
	numColumnsSubPixelMap = numColumnsPixelMap * numSubPixelsPerPixel;	// TODO Add edge pixels

	numEdgePixels = 0;
 }











/**
 * \brief: Generate the (random) flatfield variations.  This map is generated
 *		   at sub-pixel level but without the edge pixels.
 *
 * https://github.com/python-acoustics/python-acoustics/blob/master/acoustics/generator.py#L108
 */
void Detector::generateFlatfieldMap()
{

	Log.info("Detector: generating flatfield map.");

	// Random number generation

	mt19937 flatfieldGenerator(flatfieldSeed);
	normal_distribution<double> flatfieldDistribution(0.0, 1.0);

	// Double the dimensions (this is necessary because of the behaviour of the Fourier transforms)
	// (this is a bit inconvenient as we are working at sub-pixel level -> to be investigated)

	int numRows = 2 * numRowsPixelMap * numSubPixelsPerPixel;
	int numColumns = 2 * numColumnsPixelMap * numSubPixelsPerPixel;

	arma::cx_fmat evenMap = arma::cx_fmat(numRows, numColumns);

	for(unsigned int row = 0; row < numRows; row++)
	{
		for(unsigned int column = 0; column < numColumns; column++)
		{
			// Fourier space: generate white noise and include 1/f dependency
			// (Note: see https://en.wikipedia.org/wiki/Pink_noise#Generalization_to_more_than_one_dimension)

			evenMap(row, column) = flatfieldDistribution(flatfieldGenerator) / (pow(row, 2) + std::pow(column, 2) + 1);
		}
	}

	// Take the real part of the inverse Fourier transform

	evenMap = arma::ifft2(evenMap);
	arma::fmat realMap = arma::real(evenMap);

	// Cut out the appropriate part

	flatfieldMap(arma::span::all, arma::span::all) = realMap(arma::span(0, numRows / 2 - 1), arma::span(0, numColumns / 2 - 1));

	// Normalise

	float minPinkNoise = flatfieldMap.min();
	float maxPinkNoise = flatfieldMap.max();


	flatfieldMap -= minPinkNoise;
	flatfieldMap /= (maxPinkNoise - minPinkNoise); // [0, 1]
	flatfieldMap *= flatfieldNoiseAmplitude;	// [0, flatfialdNoiseAmplitude]
	flatfieldMap += (1.0 - flatfieldNoiseAmplitude);

	// Save the intra-pixel flatfield in the HDF5 file

	Log.debug("Detector: writing IRNU to HDF5");

	hdf5File.writeArray("/Flatfield", "IRNU", flatfieldMap);

	// Rebin the intra-pixel flatfield to the pixel flatfield (IRNU -> PRNU)
	// and also write this array to the HDF5 outputfile. This PRNU array is not used
	// in the remainder of the simulation.

	arma::Mat<float> prnu(numRowsPixelMap, numColumnsPixelMap, arma::fill::zeros);

	for (unsigned int row = 0; row < numRowsPixelMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsPixelMap; column++)
		{
			const unsigned int beginRow = row * numSubPixelsPerPixel;
			const unsigned int beginCol = column * numSubPixelsPerPixel;
			const unsigned int endRow = (row + 1) * numSubPixelsPerPixel - 1;
			const unsigned int endCol = (column + 1) * numSubPixelsPerPixel - 1;

			prnu(row, column) = arma::accu(flatfieldMap.submat(beginRow, beginCol, endRow, endCol))
                                / (numSubPixelsPerPixel * numSubPixelsPerPixel);
		}
	}

	// Write the result to the HDF5 output file

	Log.debug("Detector: writing PRNU to HDF5");

	hdf5File.writeArray("/Flatfield", "PRNU", prnu);
}











/**
 *\brief Generate throughput map, containing for each sub-field pixel the combined throughput efficiency
 *       of vignetting, polarisation, particulate & molecular contamination, and quantum efficiency.  Each
 *       array value is a value between 0 and 1.
 * 
 * \details Because of vignetting, the stars at the edge of the FOV look dimmer than the stars close
 *          to the optical axis. If the incoming flux before vignetting at pixel (i,j) is F(i,j), 
 *          then the flux after vignetting taken into account is F(i,j) * vignettingMap(i,j).
 *          Because of contamination (both particulate and molecular) the throughput efficiency
 *          decreases over the entire FOV by the same factor.
 *          
 * \note    The throughput map is written to the HDF5 map.
 */

void Detector::generateThroughputMap()
{
	Log.info("Detector: generating throughput map.");

	double xFPmm, yFPmm;
	double angle;

	const double refAnglePolarizationRadians = deg2rad(refAnglePolarization);		// Reference angle for the polarisation efficiency [radians]
	const double cosPolarizationEfficiency = cos(polarizationEfficiency);

	const double refAngleQuantumEfficiencyRadians = deg2rad(refAngleQuantumEfficiency);		// Reference angle for the quantum efficiency [radians]
	const double cosQuantumEfficiency = cos(quantumEfficiency);

	if (includeVignetting || includePolarization || includeQuantumEfficiency)
	{
		// Loop over all pixels in the pixel map

		for (unsigned int row = 0; row < numRowsPixelMap; row++)
		{

			for (unsigned int column = 0; column < numColumnsPixelMap;
					column++)
			{

				// Pixel coordinates (in the detector) -> focal-plane coordinates

				tie(xFPmm, yFPmm) = pixelToFocalPlaneCoordinates(row, column);

				// Angular distance [radians] of the pixel from the optical axis

				angle = camera.getGnomonicRadialDistanceFromOpticalAxis(xFPmm,
						yFPmm);

				// Vignetting

				if (includeVignetting)
					throughputMap(row, column) *= pow(cos(angle), 2);

				// Polarisation

				if (includePolarization)
					throughputMap(row, column) *= cos(
							angle / refAnglePolarizationRadians
									* cosPolarizationEfficiency);// Eq. 4-11 in PLATO-DLR-PL-RP-001

				// Quantum efficiency
				// Pixel units before: [photons]
				// Pixel units after: [electrons]

				if (includeQuantumEfficiency)
					throughputMap(row, column) *= cos(
							angle / refAngleQuantumEfficiencyRadians
									* cosQuantumEfficiency);// Eq. 4-12 in PLATO-DLR-PL-RP-001
			}
		}
	}

	// Particulate contamination

	if(includeParticulateContamination)
	{
		throughputMap *= particulateContaminationEfficiency;	// Sect. 4.2.4.3 in PLATO-DLR-PL-RP-001
	}

	// Molecular contamination

	if(includeMolecularContamination)
	{
		throughputMap *= molecularContaminationEfficiency;	// Sect. 4.2.4.4 in PLATO-DLR-PL-RP-001
	}

	// Write the result to HDF5

	Log.debug("Detector: writing throughput map to HDF5");

	hdf5File.writeArray("/Throughput", "throughputMap", throughputMap);
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

	Log.info("Detector: Integrating light for exposure " + to_string(imageNr) + " with exposure time = " + to_string(exposureTime));

	integrateLight(startTime, exposureTime);

	// Include noise effects like readout noise, photon noise, full well saturation, etc.
	// Note: readOut() needs the exposure time to compute the open shutter smearing.

	Log.info("Detector: Adding noise effects to exposure " + to_string(imageNr));

	readOut(exposureTime);

	// Write the CCD subfield, the bias map, and the smearing map to the HDF5 file

	Log.debug("Detector: Writing PixelMap, smearing map, and bias map #" + to_string(imageNr) + " to HDF5 file.");

	writePixelMapsToHDF5();

	// If required, also write the subpixel image to the HDF5 file

	if (writeSubPixelImagesToHDF5)
	{
		Log.debug("Detector: Writing SubPixelMap " + to_string(subPixelImageNr) + " to HDF5 file.");
		writeSubPixelMapToHDF5();
	}

	// Advance the internal clock

	internalTime += exposureTime + readoutTime;

	return internalTime;
}












/**
 * \brief: During an exposure, this method makes the detector integrate the light
 *         in small steps. During each step the slight change of star positions due
 *         to spacecraft jitter is taken into account. 
 *         
 *  \details  Besides jitter, also the sky background, and the flatfield is taken into 
 *            account. The sub-pixel map is rebinned in a pixel map.  After rebinning,
 *            vignetting and polarisation are applied (if applicable).
 *
 * \note The convolution with the PSF is not yet done here.
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

	Log.debug("Detector: resetting subfield array for new exposure.");

	reset();

	// Integration (incl. jitter) + background

	camera.exposeDetector(*this, startTime, exposureTime);

	// Apply flatfield (at sub-pixel level)

    if (includeFlatfield)
    {
        Log.debug("Detector: applying Flatfield.");

        applyFlatfield();
    }
    else
    {
        Log.debug("Detector: no flatfield applied.");
    }

	// Rebin from a subpixel map to a pixel map

	Log.debug("Detector: rebinning sub-pixel map into pixel map.");

	rebin();

	// Apply throughput efficiency on the pixel map

	applyThroughputEfficiency();
}












/**
 * \brief: Add the given flux value to the value of the sub-pixel that corresponds to the given coordinates 
 *         in the focal plane. Return the pixel coordinates of the pixel to which the flux was added.
 *
 * \param xFPprime   X-coordinate of the sub-pixel in the focal plane in the FP' reference frame [mm].
 * \param yFPprime   Y-coordinate of the sub-pixel in the focal plane in the FP' reference frame [mm].
 * \param flux       Flux to add to the sub-pixel map [photons].
 *
 * \return           (isInSubfield, row, col) 
 *                   isInSubfield: True if (xFPprime, yFPprime) are on the subfield, false otherwise.
 *                   row:          subfield (not CCD) row number of the pixel to which the flux was added
 *                   col:          subfield (not CCD) column number of the pixel to which the flux was added  
 */

tuple<bool, double, double> Detector::addFlux(double xFPprime, double yFPprime, double flux)
{
	// Convert from FP' coordinates to CCD pixel coordinates

	double pixRow, pixColumn;
	tie(pixRow, pixColumn) = focalPlaneToPixelCoordinates(xFPprime, yFPprime);

	// Sub-field coordinates, taking into account the edge pixels 
	// (subpixRow, subpixColumn) are the indices of the star in the subpixelMap. So they are not 
	// subpixel coordinates in the CCD frame, but in the subfield reference frame.

	const double subpixColumn = round((pixColumn - subFieldZeroPointColumn + numEdgePixels) * numSubPixelsPerPixel);
	const double subpixRow    = round((pixRow    - subFieldZeroPointRow    + numEdgePixels) * numSubPixelsPerPixel);

	// Convert back the _rounded_ subpixel coordinates to pixel coordinates
	// E.g. if there are 4 subpixels per pixel, then the pixel coordinates should always end with
	//      0.0, 0.25, 0.5, or 0.75

	pixRow    = subpixRow    / numSubPixelsPerPixel - numEdgePixels;
	pixColumn = subpixColumn / numSubPixelsPerPixel - numEdgePixels;

	// Add the flux to the subPixelMap

	if (isInSubPixelMap(subpixRow, subpixColumn))
	{
		subPixelMap((int) subpixRow, (int) subpixColumn) += flux;
		return make_tuple(true, pixRow, pixColumn);
	}
	else
	{
		return make_tuple(false, pixRow, pixColumn);
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

	double rowUnrot = (xFPprime - originOffsetY) / (pixelSize / 1000.0);
	double colUnrot = (yFPprime - originOffsetX) / (pixelSize / 1000.0);

	// Compute the coordinates in the rotated CCD reference frame

	double colRot = colUnrot * cos(orientationAngle) - rowUnrot * sin(orientationAngle);
	double rowRot = colUnrot * sin(orientationAngle) + rowUnrot * cos(orientationAngle);

	// Check wether these pixel coordinates falls on the subfield

	return    (colRot >= subFieldZeroPointColumn) && (colRot < subFieldZeroPointColumn + numColumnsPixelMap)
	       && (rowRot >= subFieldZeroPointRow)    && (rowRot < subFieldZeroPointRow + numRowsPixelMap);
}


















/**
 * \brief   Check whether the given (row, column) indices are within the array range of the subpixel map.
 *
 * \details  The input parameters row & column come from a coordinate transformation
 *           in the focal plane, and as a result are not necessarily integers. For this 
 *           function it's not necessary to round them to the nearest integer. 
 *
 * \param  row:    Row index. NOT a coordinate in the CCD frame, but in the subfield frame. [sub-pixel].
 * \param  column: Column index.NOT a coordinate in the CCD frame, but in the subfield frame.  [sub-pixel].
 *
 * \return  True if the given (row, column) coordinates are in the sub-pixel map; false otherwise.
 */

bool Detector::isInSubPixelMap(double row, double column)
{
	return (column >= 0) && (row >= 0) && (column < numColumnsSubPixelMap) && (row < numRowsSubPixelMap);
}










/**
 * \brief: Add the given flux value to (all sub-pixels of) the sub-pixel map.
 *
 * \param flux: Flux to add to the sub-pixel map [photons/pixel].
 *
 */

void Detector::addFlux(double flux)
{
	// The flux is expressed in [photons/pixel] but we need the quantity expressed 
	// in [photons/subpixel]. There are (numSubPixelsPerPixel)^2 per pixel (the
	// name is thus a bit of a misnomer.).

	subPixelMap += flux / numSubPixelsPerPixel / numSubPixelsPerPixel;
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
 * \brief Apply throughput efficiency. This is the combined effect of:
 * 			- vignetting (brightness attenuation towards the edges of the FOV);
 * 			- particulate contamination;
 * 			- molecular contamination;
 * 			- quantum efficiency.
 */

void Detector::applyThroughputEfficiency()
{
    pixelMap = pixelMap % throughputMap;		// Element-wise multiplication with the throughput map
}














/**
 * \brief: Reads out the detector and apply the following effects:
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

	// Apply poisson distributed photon noise
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	if (includePhotonNoise)
	{
        Log.debug("Detector: adding photon noise");
		addPhotonNoise();
	}
    else 
    {
        Log.debug("Detector: no photon noise added.");
    }

	// Apply full-well saturation. A pixel has a maximum capacity of electrons (the full well capacity).
	// If photons free more electrons, the pixel saturates, and the electrons flow in the pixels above and below in
	// the same column (potential barriers are smallest in that direction).
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

    if (includeFullWellSaturation)
    {
        Log.debug("Detector: applying full well saturation.");
        applyFullWellSaturation();
    }
    else
    {
        Log.debug("Detector: no full well saturation applied.");
    }

	// Simulate the effects of the Charge Transfer Inefficiency (CTI). When the
	// CCD is read out, row after row, a part of the charge is always left behind
	// which then dribbles into the trailing pixels. This causes each star to have
	// a small "tail". Only visible when the CTI = 1 - CTE is poor.
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	if (includeCTIeffects)
	{
        Log.debug("Detector: applying charge transfer inefficiency");
		applyCTI();
	}
    else
    {
        Log.debug("Detector: no charge transfer inefficiency applied.");
    }

	// Apply the effects of readout smearing due to an open shutter. Because there is no shutter,
	// the pixels are still receiving photons from the sky, while they are being transfered towards
	// the readout register.

	if (includeOpenShutterSmearing)
	{
        Log.debug("Detector: applying open shutter smearing.");
		applyOpenShutterSmearing(exposureTime);
	}
    else 
    {
        Log.debug("Detector: no open shutter smearing applied.");
    }

	// Each time the amplifier reads out a pixel, a tiny bit of noise is added.
	// Add the readout noise.
	// Pixel units before: [electrons]
	// Pixel units after: [electrons]

	if (includeReadoutNoise)
	{ 
        Log.debug("Detector: adding readout noise.");
		addReadoutNoise();
	}
    else
    {
        Log.debug("Detector: no readout noise added.");
    }

	// Apply the gain, to increase the dynamic range of the detector.
	// Pixel units before: [electrons]
	// Pixel units after: [ADU]

	applyGain();

	// Take into account the bias level (i.e. add the constant "zero" level
	// introduced by the amplifier).
	// Pixel units before: [ADU]
	// Pixel units after: [ADU]

	addElectronicOffset();


    // At this point, because pixel values are floats, the previous steps may have
    // resulted in fractional ADUs. Take care that pixel maps, bias maps, and smearing
    // maps do not have fractional values.

    pixelMap = arma::floor(pixelMap);
    biasMap = arma::floor(biasMap);
    smearingMap = arma::floor(smearingMap);
    

	// Take into account digital saturation. If even after dividing by the gain
	// the number of ADUs in a pixel is still higher than the analogue-digital
	// converter (ADC) can represent with its fixed amount of bits, clip all
	// values that are too high to the saturation level of the ADC.
	// Pixel units before: [ADU]
	// Pixel units after: [ADU]

    if (includeDigitalSaturation)
    { 
        Log.debug("Detector: applying digital saturation to pixelMap, biasMap and smearingMap (digitalSaturationLimit=" + to_string(digitalSaturationLimit) + ")");
    	applyDigitalSaturation();
    }
    else
    {
        Log.debug("Detector: no digital saturation applied.");
    }
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
				// Transfer excess electrons up

				jmod = row;
				numExcessElectrons = (pixelValue - fullWellSaturationLimit) / 2.0;   // Move half of the excess electrons down...

                // When we move half of the excess electrons down, the pixel below may also become saturated 
                // (or was already saturated). When processing this pixel below, we need to avoid that it also 
                // sends its excess electrons up and down, and therefore sends part of its excess electrons 
                // back up where they actually came from. The variable 'transfer2Saturated' is to keep track 
                // whether the destination pixel is also saturated, so that it does not send back electrons
                // where they came from.
                
                bool transfer2Saturated = false;

				while (numExcessElectrons > 0 && jmod < numRowsPixelMap)
				{
					if(!transfer2Saturated)
					{
						pixelMap(jmod, column) -= numExcessElectrons;
					}

					jmod++;

					// Electrons reaching the edge of the CCD will not be detected

					if (jmod < numRowsPixelMap)
					{
						if(pixelMap(jmod, column) >= fullWellSaturationLimit)
						{
							transfer2Saturated= true;
						}

						else{

							transfer2Saturated = false;

							pixelMap(jmod, column) += numExcessElectrons;

							// Make sure the pixel you move the excess electrons to
							// does not get saturated too

							if (pixelMap(jmod, column) > fullWellSaturationLimit)
							{
								numExcessElectrons = pixelMap(jmod, column) - fullWellSaturationLimit;
							}

							else
							{
								numExcessElectrons = 0;
							}
						}
					}
				}

				// Transfer excess electrons down

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

						else
						{
							numExcessElectrons = 0;
						}
					}
				}
			}
		}
	}
}











/**
 * \brief Apply the effect of the charge-transfer inefficiency to the
 *        pixel map. The exact model used depends on the configuration
 *        in the input file.
 *  
 *  \note The pixel map should be expressed in [e-] and not [ADU]
 *  
 * \pre Pixel unit in the pixel map: [electrons]
 * \pre Pixel unit in the smearing map: [electrons].
 * \pre No bias register map.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map.              
 */

void Detector::applyCTI()
{
    if (CTImodel == "Simple")
    {
        Log.info("Detector: applying simple CTI model");
        applySimpleCTImodel();
    }
    else if (CTImodel == "Short2013")
    {
        Log.info("Detector: applying Short et al. (2013) CTI model");
        applyShort2013CTImodel();
    }
}











/**
 * \brief: Apply the effect of the charge-transfer inefficiency to the
 *         pixel map, using a simple CTI model. The CTI of this simple model
 *         has no dependence on the flux level, nor on the distance of the 
 *         readout register.
 *         
 * \note The serial register is assumed to have a CTE of 1.
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre Pixel unit in the smearing map: [electrons].
 * \pre No bias register map.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map.
 */

void Detector::applySimpleCTImodel()
{
	float cti = 1.0 - meanCte;

	// Computing the effects of CTE requires the use of a binomial distribution.
	// To speed up things, we first pre-compute some parts of this distribution.

	// Pre-compute the (natural) logarithms of the first N natural numbers

	vector<double> logs(numRowsPixelMap + subFieldZeroPointRow);
	iota(logs.begin(), logs.end(), 1.0);
	transform(logs.begin(), logs.end(), logs.begin(), ptr_fun<double, double>(log));

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

		if (row + subFieldZeroPointRow == 0)
		{
			const double factor1 = meanCte;

			readout += pixelMap(0, arma::span::all) * factor1;

		}
		else
		{
			for (unsigned int index = subFieldZeroPointRow; index <= row + subFieldZeroPointRow; index++)
			{
				const double cteFactor = pow(meanCte, index + 1)
						* pow(cti, row + subFieldZeroPointRow - index);

				if ((index == 0) || (row - (index - subFieldZeroPointRow) == 0))
				{
					readout += pixelMap(index - subFieldZeroPointRow, arma::span::all) * cteFactor;
				}

				else
				{
					const double binomialFactor = exp(sumOfLogsUpTo[row + subFieldZeroPointRow - 1]
									                  - sumOfLogsUpTo[row - (index - subFieldZeroPointRow) - 1]
									                  - sumOfLogsUpTo[index - 1]);

					readout += pixelMap(index - subFieldZeroPointRow, arma::span::all) * cteFactor;
				}
			}
		}

		pixelMap(row, arma::span::all) = readout(0, arma::span::all);
	}
}











/**
 * \brief: Apply the effect of the charge-transfer inefficiency to the pixel map,
 *         using the model described in Short et al., MNRAS 430, 3078-3085 (2013).
 *         Only parallel readout is taken into account here.
 *         
 * \note The readout register is assumed to be right next to row [0] of the pixel map.
 *       The pixel map needs to be in [e-], not in [ADU]
 * 
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre Pixel unit in the smearing map: [electrons].
 * \pre No bias register map.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map.
 */
 
void Detector::applyShort2013CTImodel()
{
    // Compute the maximum geometrical volume that electrons can occupy within a pixel.
    // I.e. the volume of the electron cloud when the capacity of the full well is maxed out.
    // E.g if the pixel size is 18 micron, then the volume is 18e-6 * 18e-6 * 1.e-6 / 2.0

    const double maxVolumePerPixel = pixelSize * pixelSize * 1.e-18  / 2.0;                                  // [m^3]

    // Compute the time it takes to transfer 1 row during readout

    const double chargeTransferTime = readoutTime / numRows;                                                 // [s] 

    // Compute the thermal velocity of the electrons in the silicon

    const double effectiveElectronMass = 0.5 * Constants::FREEELECTRONMASS;                                  // [kg]
    const double thermalVelocity = sqrt(3.0 * Constants::KBOLTZMANN * temperature / effectiveElectronMass);  // [m/s]

    // Arrays to keep track of the number of occupied traps in a column

    arma::Mat<float> numberOfOccupiedTraps = arma::zeros<arma::Mat<float>>(numTrapSpecies, numColumnsPixelMap);

    // Arrays to keep track of the captured and released electrons, for each column in a particular row.

    arma::Row<float> numberOfCapturedElectrons(numColumnsPixelMap);
    arma::Row<float> numberOfReleasedElectrons(numColumnsPixelMap);

    // Loop over all rows of the pixelMap, and over all trap species.
    // For each row, the computations are done for all columns simultaneously.

    for (int rowNumber = 0; rowNumber < numRowsPixelMap; rowNumber++)
    {
        for (int k = 0; k < numTrapSpecies; k++)
        {
            // Compute the number of electrons captured in a trap, according to Eq. (22)-(23) of Short et al. (2013).
            // Note that Armadillo uses % for elementwise multiplication.

            const double alpha = chargeTransferTime * trapCaptureCrossSection[k] * thermalVelocity * pow(fullWellSaturationLimit, beta) / 2.0 / maxVolumePerPixel;
            const double gamma = trapDensity[k] * (subFieldZeroPointRow + rowNumber) / pow(fullWellSaturationLimit, beta);

            numberOfCapturedElectrons =   (gamma * arma::pow(pixelMap.row(rowNumber), beta) - numberOfOccupiedTraps.row(k)) \
                                        / (gamma * arma::pow(pixelMap.row(rowNumber), beta-1) + 1)                          \
                                        % (1 - arma::exp(-alpha * arma::pow(pixelMap.row(rowNumber), 1-beta)));

            // Captured electron numbers can't be negative, so clip negative value to zero.

            arma::Col<arma::uword> isNegative = arma::find(numberOfCapturedElectrons < 0.0);
            numberOfCapturedElectrons(isNegative).zeros();

            // Update the number of occupied traps with the estimated number of captured electrons

            numberOfOccupiedTraps.row(k) += numberOfCapturedElectrons;
            
            // Correct the number of occupied traps with the electrons that were released again during the charge transfer time.

            numberOfReleasedElectrons = numberOfOccupiedTraps.row(k) * (1-exp(-chargeTransferTime/releaseTime[k]));
            numberOfOccupiedTraps.row(k) -= numberOfReleasedElectrons;

            // Add the electron excess to the current pixel value

            pixelMap.row(rowNumber) += numberOfReleasedElectrons - numberOfCapturedElectrons;
        }
    }
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
	// Average out the fluxes in the pixel map per column (of the whole CCD) and make sure it is
	// scaled with the readout time instead of with the exposure time:
	// - rows in the sub-field: use actual fluxes
	// - rows outside the sub-field: use total sky background

	arma::Row<float> openShutterSmearing = arma::sum(pixelMap, 0);

	float openShutterSmearingOutsideSubField = camera.getTotalSkyBackground() * (numRows - numRowsPixelMap + numRowsBiasMap + numRowsSmearingMap);

	// Apply all throughput efficiencies

	if(includeVignetting)
		openShutterSmearingOutsideSubField *= expectedValueVignetting;

	if(includePolarization)
		openShutterSmearingOutsideSubField *= expectedValuePolarization;

	if(includeParticulateContamination)
		openShutterSmearingOutsideSubField *= particulateContaminationEfficiency;

	if(includeMolecularContamination)
		openShutterSmearingOutsideSubField *= molecularContaminationEfficiency;

	if(includeQuantumEfficiency)
		openShutterSmearingOutsideSubField *= expectedValueQuantumEfficiency;


	openShutterSmearing += openShutterSmearingOutsideSubField;

	float factor = (readoutTime / exposureTime) / numRows;
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
 * \brief Apply the readout noise to the pixel map, bias map, and smearing map
 * 
 * \details Readout noise occurs due to the imperfect nature of the CCD amplifiers.  
 *          When the electrons are transferred to the amplifier, the induced voltage
 *          is measured. However, this measurement is not perfect, but gives a value 
 *          which is on average correct, with the readout noise as standard deviation.
 *          So readout noise is a measure of this scatter around the true value.
 *          Its value is expressed in electrons as the packet of charge is made up of 
 *          electrons.
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre Pixel unit in the smearing map: [electrons].
 * \pre Bias map not initialised.
 * \pre Readout noise expressed in [electrons].
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post Pixel unit in the bias register map: [electrons].
 * \post Initialised the bias and smearing maps with readout noise.
 */

void Detector::addReadoutNoise()
{

	readoutNoiseDistribution = normal_distribution<double>(0.0, readoutNoise);

	// Add readout noise to the pixel map

	for (unsigned int row = 0; row < numRowsPixelMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsPixelMap; column++)
		{
			pixelMap(row, column) += readoutNoiseDistribution(readoutNoiseGenerator);
		}
	}

	// Add readout noise to the bias prescan map

	for (unsigned int row = 0; row < numRowsBiasMap; row++)
	{
		for (unsigned int column = 0; column < numColumnsPixelMap; column++)
		{
			biasMap(row, column) += readoutNoiseDistribution(readoutNoiseGenerator);
		}
	}

    // Add readout noise to the smearing overscan map

    for (unsigned int row = 0; row < numRowsSmearingMap; row++)
    {
        for (unsigned int column = 0; column < numColumnsPixelMap; column++)
        {
            smearingMap(row, column) += readoutNoiseDistribution(readoutNoiseGenerator);
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
	Log.debug("Detector: applying gain to pixelMap, biasMap and smearingMap (gain=" + to_string(gain) + ")");

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
	Log.debug("Detector: adding a bias to pixelMap, biasMap and smearingMap (electronicOffset=" + to_string(electronicOffset)+ ")");

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
	// Top off the values in the pixel map

	pixelMap(arma::find(pixelMap > digitalSaturationLimit)).fill(digitalSaturationLimit);

	// Top off the values in the bias register map

	biasMap(arma::find(biasMap > digitalSaturationLimit)).fill(digitalSaturationLimit);

	// Top off the values in the smearing map

	smearingMap(arma::find(smearingMap > digitalSaturationLimit)).fill(digitalSaturationLimit);
}












/**
 * \brief Compute the (x,y) coordinates [mm] in the FP' reference system (not the FP system) 
 *        given the (real-valued) pixel row and column numbers on the CCD.
 *        
 * \note  The rows correspond to the y-direction, and the columns to the x-direction.
 *        Pixel (row, col) = (0,0) starts at (yFP, xFP) = (0, 0).
 *               
 * \param row     Row coordinate, real-valued (e.g. 3.5)    [pix]
 * \param column  Column coordinate, real-valued (e.g. 8.3) [pix]
 * 
 * \return (xFPprime, yFPprime)  A pair of (x,y) coordinates in the FP' reference system [mm]
 */

pair<double, double> Detector::pixelToFocalPlaneCoordinates(double row, double column)
{
    // Convert the pixel coordinates into [mm] coordinates
    // The pixelSize is expressed in [micron].

    double xCCDmm = column * pixelSize / 1000.0;
    double yCCDmm = row * pixelSize / 1000.0;

    // Convert the CCD coordinates into FP' coordinates [mm]
    // Note: orientationAngle is in [rad], originOffsetX and originOffsetY in mm

    double xFPprime = originOffsetX + xCCDmm * cos(orientationAngle) - yCCDmm * sin(orientationAngle);
    double yFPprime = originOffsetY + xCCDmm * sin(orientationAngle) + yCCDmm * cos(orientationAngle);

    // That's it

    return make_pair(xFPprime, yFPprime);
}










/**
 * \brief Compute the (real-valued) pixel coordinates of the star on the CCD, given the 
 *        (x,y) coordinates [mm] in the FP' reference system (not the FP system)
 *
 * \note  The rows correspond to the y-direction, and the columns to the x-direction.
 *        Pixel (row, col) = (0,0) starts at (yFP, xFP) = (0, 0).
 *        
 * \param xFPprime  x-coordinate of the point in the FP' reference system  [mm]
 * \param yFPprime  y-coordinate of the point in the FP' reference system  [mm]
 * 
 * \return (row, column)  Row and column pixel coordinates of the point (real-valued) [pix]
 */

pair<double, double> Detector::focalPlaneToPixelCoordinates(double xFPprime, double yFPprime)
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
 * \return (xFPprime, yFPprime)   focal plane coordinates in the FP' reference system [mm]
 */

pair<double, double> Detector::getFocalPlaneCoordinatesOfSubfieldCenter()
{
	double centerRow = subFieldZeroPointRow + numRowsPixelMap / 2.0;
	double centerCol = subFieldZeroPointColumn + numColumnsPixelMap / 2.0;

    // The columns correspond to the x-coordinate, the rows to the y-coordinate

    double xFPprime, yFPprime;
    tie(xFPprime, yFPprime) = pixelToFocalPlaneCoordinates(centerRow, centerCol);

	return make_pair(xFPprime, yFPprime);
}












/**
 * \brief Return a boolean telling whether the PSF has been previously set, e.g.
 *        through setPsfForSubfield().
 * 
 * \return  true if the PSF was previously set, false otherwise
 */

bool Detector::psfIsSet()
{
	return psfWasSet;
}









/**
 * \brief      Set the Point Spread Function map for the subfield
 * 
 * \details    The PSF that is selected is dependent on the user input.
 */
void Detector::setPsfForSubfield()
{
    double centerXmm, centerYmm;
    tie(centerXmm, centerYmm) = getFocalPlaneCoordinatesOfSubfieldCenter();

    arma::fmat psf = camera.getRebinnedPsfForFocalPlaneCoordinates(centerXmm, centerYmm, numSubPixelsPerPixel, getOrientationAngle());

	convolver.initialise(numRowsSubPixelMap, numColumnsSubPixelMap, psf);

    psfMap = psf;
    psfWasSet = true;
}











/**
 * \brief: Convolve the sub-pixel map with the PSF, keeping the same dimensions.
 *
 * \param psf: PSF.
 */
void Detector::convolveWithPsf()
{

    if(includeConvolution)
    {
        Log.debug("Detector: convolving subPixelMap with PSF.");

        // subpixelMap serves here both as input as well as output matrix;

    	convolver.convolve(subPixelMap, subPixelMap);
    }
    else
    {
        Log.debug("Detector: no convolution applied.");
    }

}













/**
 * \brief Return the (X,Y) coordinates in the FP' reference frame in [mm] of the 4 corners
 *        of the subfield
 *        
 * \return (X00, Y00, X01, Y01, X11, Y11, X10, Y10)  [mm]
 *         where: (X00, Y00) are the FP' coordinates of the lower left corner of the subfield
 *                (X01, Y01) are the FP' coordinates of the lower right corner of the subfield
 *                (X11, Y11) are the FP' coordinates of the upper right corner of the subfield
 *                (X10, Y10) are the FP' coordinates of the upper left corner of the subfield
 */

tuple<double, double, double, double, double, double, double, double> Detector::getFocalPlaneCoordinatesOfSubfieldCorners()
{
	double corner00Xmm, corner00Ymm, corner01Xmm, corner01Ymm, corner11Xmm, corner11Ymm, corner10Xmm, corner10Ymm;
	double row, col;

	// Lower left corner

	row = subFieldZeroPointRow;
	col = subFieldZeroPointColumn;
	tie(corner00Xmm, corner00Ymm) = pixelToFocalPlaneCoordinates(row, col);

	// Lower right corner

	row = subFieldZeroPointRow;
	col = subFieldZeroPointColumn + numColumnsPixelMap;
	tie(corner01Xmm, corner01Ymm) = pixelToFocalPlaneCoordinates(row, col);

	// Upper right corner

	row = subFieldZeroPointRow + numRowsPixelMap;
	col = subFieldZeroPointColumn + numColumnsPixelMap;
	tie(corner11Xmm, corner11Ymm) = pixelToFocalPlaneCoordinates(row, col);

	// Upper left corner

	row = subFieldZeroPointRow + numRowsPixelMap;
	col = subFieldZeroPointColumn;
	tie(corner10Xmm, corner10Ymm) = pixelToFocalPlaneCoordinates(row, col);

	return make_tuple(corner00Xmm, corner00Ymm, corner01Xmm, corner01Ymm, corner11Xmm, corner11Ymm, corner10Xmm, corner10Ymm);
}












/**
 * \brief Return the solid angle of 1 single pixel on the sky. [sr]
 * 
 * \param plateScale  The platescale of the camera [arcsec/micron]
 *
 * \return            Solid angle in [s]
 */

double Detector::getSolidAngleOfOnePixel(double plateScale)
{
	return sqDeg2sr(pow(pixelSize * plateScale / 3600.0, 2));
}












/**
 * \brief      Return the orientation of the CCD with respect to the orientation
 *             of the focal plane. The rotations of the CCD are counter
 *             clockwise.
 *
 * \return     the orientation of the CCD [radians]
 */
double Detector::getOrientationAngle()
{
    return orientationAngle;
}












/**
 * \brief     Set the subfield with a given array.  
 * 
 * \details   This function is primarily used for testing the code. One should not first get the pixelMap
 *            perform an operation, and then setSubfield() again. Instead, let Detector do the operation.
 *  
 * \param subfield
 */

void Detector::setSubfield(const arma::Mat<float> &subfield)
{
	// Check if the given matrix has the proper dimensions. If not complain, and exit.

	if ((subfield.n_rows != pixelMap.n_rows) || (subfield.n_cols != pixelMap.n_cols))
	{
		Log.error("Detector: setSubfield with incompatible array shape: (" 
		          + to_string(subfield.n_rows) + ", " + to_string(subfield.n_cols) + ") != ("
		          + to_string(pixelMap.n_rows) + ", " + to_string(pixelMap.n_cols) + ")");
		exit(1);
	} 

	// Copy the contents of the subfield array into our pixelMap

	pixelMap = subfield;
}












/**
 * \brief Return a copy of the pixelMap matrix
 * 
 * \details   This function is primarily used for testing the code. One should not first get the pixelMap
 *            perform an operation, and then setSubfield() again. Instead, let Detector do the operation.
 * 
 * \return pixelMap
 */

 arma::Mat<float> Detector::getSubfield()
 {
 	return pixelMap;
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
	hdf5File.createGroup("/Throughput");

	if (writeSubPixelImagesToHDF5)
	{
		hdf5File.createGroup("/SubPixelImages");
	}
}












/**
 * \brief: Writes the pixel map for the HDF5 file.
 */

void Detector::writePixelMapsToHDF5()
{
    // Compose the image name

	stringstream myStream;
    myStream << "image" << setfill('0') << setw(6) << imageNr;
    string imageName = myStream.str();

    // Add the image to the "Images" group

    hdf5File.writeArray("/Images", imageName, pixelMap);

    // Clear the string stream and compose the smearing map name

    myStream.str(string());      // insert empty string
    myStream.clear();            // clear eof bit

    myStream << "smearingMap" << setfill('0') << setw(6) << imageNr;
    string smearingMapName = myStream.str();

    // Add the smearing map to the "SmearingMaps" group

    hdf5File.writeArray("/SmearingMaps", smearingMapName, smearingMap);

    // Clear the string stream and compose the bias map name

    myStream.str(string());      // insert empty string
    myStream.clear();            // clear eof bit

    myStream << "biasMap" << setfill('0') << setw(6) << imageNr;
    string biasMapName = myStream.str();

    // Add the smearing map to the "SmearingMaps" group

    hdf5File.writeArray("/BiasMaps", biasMapName, biasMap);


    // Increment the counter for the next image

    imageNr++;
}











/**
 * \brief: Writes the subpixel map for the HDF5 file.
 */

void Detector::writeSubPixelMapToHDF5()
{
	stringstream myStream;
    myStream << "subPixelImage" << setfill('0') << setw(6) << subPixelImageNr;
    string imageName = myStream.str();

    // Add the image to the "SubPixelImages" group

    hdf5File.writeArray("/SubPixelImages", imageName, subPixelMap);

    // Increment the counter for the next subpixel image

    subPixelImageNr++;
}


