#include "DetectorWithMappedPSF.h"

/**
 * \brief Default constructor.
 * 
 * \param configParam: Configuration parameters.
 * 
 * \param hdf5file: HDF5 file from which to read the PSFs.
 */

DetectorWithMappedPSF::DetectorWithMappedPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
: Detector(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure), includeFlatfield(true), writeSubPixelImagesToHDF5(false)
{
    // Parse the parameters from the configuration file.

    configure(configParam);

    // Create the groups in the HDF5 file where the different maps (i.e. pixel map,
    // bias register map, smearing map, etc.) will be saved. This needs to be done
    // BEFORE other methods write arrays to HDF5.

    initHDF5Groups();

    // Allocate memory for the different maps

    subPixelMap.zeros(numRowsSubPixelMap, numColumnsSubPixelMap);

    flatfieldMap.ones(numsubsubfieldsx * (numRowsSubPixelMap - 2*overlapx * numSubPixelsPerPixel) + 2*overlapx * numSubPixelsPerPixel, numsubsubfieldsy * (numColumnsSubPixelMap - 2*overlapy * numSubPixelsPerPixel) + 2*overlapy * numSubPixelsPerPixel);  //%% For spectral dependency, make flatfield map cover the full final stitched field

    if(includeFlatfield)
    {
        // Generate the flatfield map

        generateFlatfieldMap();
    }

    if(includeBFE)
    {
    		// Generate Guyonnet coefficients

    		generateGuyonnetCoefficients();
    }

    // Initialize and load the PSF. This will open the PSF HDF5 file and perform some basic checking, 
    // Then select the proper PSF for the given subfield. Should only be done after calling configure().

    psf = new PointSpreadFunction(configParam, hdf5file);

    for (int subsubfieldx = 0; subsubfieldx < numsubsubfieldsx; subsubfieldx++)  //%% Loop to select a different psf for each subsubfield
    {
        for (int subsubfieldy = 0; subsubfieldy < numsubsubfieldsy; subsubfieldy++)
        {
    setPsfForSubfield(subsubfieldx, subsubfieldy);  //%% Added subsubfield
        }
    }

}








/**
 * Destructor.
 *
 */
DetectorWithMappedPSF::~DetectorWithMappedPSF()
{
    flushOutput();
    delete psf;
}










/**
 * \brief Configure the DetectorWithMappedPSF object using the ConfigurationParameters
 * 
 * \param configParam: the configuration parameters 
 **/

 void DetectorWithMappedPSF::configure(ConfigurationParameters &configParam)
 {

	flatfieldNoiseRMS = configParam.getDouble("CCD/FlatfieldNoiseRMS");
	includeFlatfield = configParam.getBoolean("CCD/IncludeFlatfield");
	includeConvolution = configParam.getBoolean("CCD/IncludeConvolution");

        wave_bins = configParam.getInteger("Camera/WavelengthBins");  //%% Read how many wavelength bins are to be processed, for spectral dependency

	writeSubPixelImagesToHDF5 = configParam.getBoolean(
			"ControlHDF5Content/WriteSubPixelImages");

	numSubPixelsPerPixel = configParam.getInteger("SubField/SubPixels");

	// Configuration parameters for the noise source random seeds

	flatfieldSeed = configParam.getLong("RandomSeeds/FlatFieldSeed");

	// Treat the specific configurations for a Mapped PSF

	string psfModel = configParam.getString("PSF/Model");

	if((psfModel == "MappedGaussian") || (psfModel == "MappedFromFile"))
	{
		includeChargeDiffusion = configParam.getBoolean("PSF/" + psfModel + "/IncludeChargeDiffusion");
		includeJitterSmoothing = configParam.getBoolean("PSF/" + psfModel + "/IncludeJitterSmoothing");
		chargeDiffusionStrength = configParam.getDouble("PSF/" + psfModel + "/ChargeDiffusionStrength");

		if(includeChargeDiffusion)
		{
			generateDiffusionKernel(chargeDiffusionStrength * numSubPixelsPerPixel);
		}
		else if(includeJitterSmoothing)
		{
			generateDiffusionKernel(0.5);
		}
	}

	// Derive the dimensions of the sub-pixel map

	numRowsSubPixelMap = numRowsPixelMap * numSubPixelsPerPixel; // TODO Add edge pixels
	numColumnsSubPixelMap = numColumnsPixelMap * numSubPixelsPerPixel; // TODO Add edge pixels

}









  * \brief: Generate the diffusion kernel.  This is generated at sub-pixel level.
  *
  * \param kernelWidth: Width (sigma) of the Gaussian diffusion kernel [sub-pixels].
  */
void DetectorWithMappedPSF::generateDiffusionKernel(double kernelWidth)
{
    diffusionKernelWidth = kernelWidth;
    diffusionKernelImageSize = 2 * (int)(8 * kernelWidth + 1) + 1;

    diffusionKernel.zeros(diffusionKernelImageSize, diffusionKernelImageSize);
}

/**
 * \brief: Generate the (random) flatfield variations.  This map is generated
 *         at sub-pixel level but without the edge pixels.
 *
 * https://github.com/python-acoustics/python-acoustics/blob/master/acoustics/generator.py#L108
 */

void DetectorWithMappedPSF::generateFlatfieldMap()
{
    Log.info("DetectorWithMappedPSF: generating flatfield map.");

    // Random number generation

    mt19937 flatfieldGenerator(flatfieldSeed);
    normal_distribution<double> flatfieldDistribution(0.0, 1.0);

    // Double the dimensions (this is necessary because of the behaviour of the Fourier transforms)
    // (this is a bit inconvenient as we are working at sub-pixel level -> to be investigated)

    int Nrows = 2 * ((numsubsubfieldsx * (numRowsPixelMap - 2*overlapx) + 2*overlapx) * numSubPixelsPerPixel);  //%% For spectral dependency, FF map is size of large stutched map
    int Ncolumns = 2 * ((numsubsubfieldsy * (numColumnsPixelMap - 2*overlapy) + 2*overlapy) * numSubPixelsPerPixel);  //%% For spectral dependency, FF map is size of large stutched map
  
    arma::cx_fmat evenMap = arma::cx_fmat(Nrows, Ncolumns);

    for (unsigned int row = 0; row < Nrows; row++)
    {
        for (unsigned int column = 0; column < Ncolumns; column++)
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

    unsigned int numRowsFlatfield = Nrows / 2;
    unsigned int numColumnsFlatfield = Ncolumns / 2;

    flatfieldMap(arma::span::all, arma::span::all) = realMap(arma::span(0, numRowsFlatfield - 1), arma::span(0, numColumnsFlatfield - 1));
    flatfieldMap.reshape(numRowsFlatfield * numColumnsFlatfield, 1);

    // Normalisation
    //  - divide by mean and subtract 1.0 -> mean = 0.0
    //  - scale such that std.dev. = flatfield RMS and mean = 0.0
    //  - add 1.0

    flatfieldMap /= arma::mean(flatfieldMap.col(0));
    flatfieldMap -= 1;
    double scale = flatfieldNoiseRMS / arma::stddev(flatfieldMap.col(0));
    flatfieldMap *= scale;
    flatfieldMap += 1;

    flatfieldMap.reshape(numRowsFlatfield, numColumnsFlatfield);

    // Write the result to the HDF5 output file

    if (writeFlatfieldMap)
    {
        Log.debug("DetectorWithMappedPSF: writing IRNU to HDF5");
        hdf5File.writeArray("/Flatfield", "IRNU", flatfieldMap);
    }

    // Rebin the intra-pixel flatfield to the pixel flatfield (IRNU -> PRNU)
    // and also write this array to the HDF5 outputfile. This PRNU array is not used
    // in the remainder of the simulation.

    arma::Mat<float> prnu(numRowsPixelMap2, numColumnsPixelMap2, arma::fill::zeros);  //%% Changed to pixelMap2 (large) for spectral dependency

    for (unsigned int row = 0; row < numRowsPixelMap2; row++)  //%%
    {
        for (unsigned int column = 0; column < numColumnsPixelMap2; column++)  //%%
        {
            const unsigned int beginRow = row * numSubPixelsPerPixel;
            const unsigned int beginCol = column * numSubPixelsPerPixel;
            const unsigned int endRow = (row + 1) * numSubPixelsPerPixel - 1;
            const unsigned int endCol = (column + 1) * numSubPixelsPerPixel - 1;

            prnu(row, column) = arma::accu(flatfieldMap.submat(beginRow, beginCol, endRow, endCol)) / (numSubPixelsPerPixel * numSubPixelsPerPixel);
        }
    }

    // Write the result to the HDF5 output file

    if (writeFlatfieldMap)
    {
        Log.debug("DetectorWithMappedPSF: writing PRNU to HDF5");
        hdf5File.writeArray("/Flatfield", "PRNU", prnu);
    }
}

/**
 * \brief: Zeroes the pixel, bias register, and the smearing maps.
 *
 * \pre pixel, bias register, and smearing maps filled with values from previous exposure.
 *
 * \post pixel, bias register, and smearing maps filled with zeroes.
 */

void DetectorWithMappedPSF::reset()
{
    pixelMap.zeros();
    biasMapLeft.zeros();
    biasMapRight.zeros();
    smearingMap.zeros();
    subPixelMap.zeros();
}

/**
 * \brief: Take an exposure with the detector starting at the given time.
 *         The light is integrated during the given exposure time, during which 
 *         the detector experiences the effects of jitter and thermo-elastic telescope 
 *         drift. The background is assumed uniform for the whole subfield.
 *         Afterwards, the collected light is read out, convolving the image with the
 *         point spread function and adding various noise effects.
 *
 * \param exposureNr: Sequential number of the exposure.
 * 
 * \param startTime: Starting time of the exposure [s].
 * 
 * \param exposureTime: Duration of the exposure [s].
 * 
 * \return endTime: Time after the exposure (startTime + exposureTime + readoutTime)
 *
 * \pre Sub-pixel, pixel, bias register, and smearing map filled with values from previous exposure.
 *
 * \post Pixel unit in the pixel, bias register, and smearing maps: [ADU]
 */
double DetectorWithMappedPSF::takeExposure(int exposureNr, double startTime, double exposureTime)
{
    pixelMap2.zeros(); 	//%% Added for spectral dependency, larger map to add all subsubfields to

    // Advance the internal clock until the given start time

    internalTime = startTime;

    // Integration of point sources and background, taking into account jitter + drift.

    Log.info("DetectorWithMappedPSF: Integrating light for exposure " + to_string(exposureNr) + " with exposure time = " + to_string(exposureTime));

    for (int subsubfieldx = 0; subsubfieldx < numsubsubfieldsx; subsubfieldx++)  //%% For spectral dependency: loop over all subfields and bins
    {
        for (int subsubfieldy = 0; subsubfieldy < numsubsubfieldsy; subsubfieldy++)
        {
	integrateLight(exposureNr, startTime, exposureTime, subsubfieldx, subsubfieldy);
        }
    }

    if(includeDarkSignal)	//%% Dark signal only applies once, not wavelength dependent
    {
	Log.debug("Detector: adding dark current");

       	addDarkSignal(exposureTime);
    }
    else
    {
        Log.debug("Detector: no dark current added");
    }

    // Include noise effects like readout noise, photon noise, full well saturation, etc.
    // Note: readOut() needs the exposure time to compute the open shutter smearing.

    Log.info("DetectorWithMappedPSF: Adding noise effects to exposure " + to_string(exposureNr));

    readOut(exposureTime);

    // Write the CCD subfield, the bias map, and the smearing map to the HDF5 file

    Log.debug("DetectorWithMappedPSF: Writing PixelMap, smearing map, and bias map #" + to_string(exposureNr) + " to HDF5 file.");

    writePixelMapsToHDF5(exposureNr);

    // If required, also write the subpixel image to the HDF5 file

    if (writeSubPixelImagesToHDF5)
    {
        Log.debug("DetectorWithMappedPSF: Writing SubPixelMap " + to_string(exposureNr) + " to HDF5 file.");
        writeSubPixelMapToHDF5(exposureNr);
    }

    // Advance the internal clock

    internalTime += exposureTime + readoutTimeBeforeNextExposure;

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
 * \param exposureNr: Sequential number of the exposure.
 * 
 * \param startTime: Starting time of the exposure for which jitter must be applied [s].
 * 
 * \param exposureTime: Duration of the exposure [s].
 *
 * \pre Sub-pixel, pixel, bias register, and smearing map filled with values from previous exposure.
 *
 * \post Pixel unit of the sub-pixel map: [photons].
 * 
 * \post Pixel, bias register, and smearing map filled with zeroes.
 */

void DetectorWithMappedPSF::integrateLight(int exposureNr, double startTime, double exposureTime, int subsubfieldx, int subsubfieldy)
{
    // Reset the sub-field (i.e. get rid of the previous exposure, by zeroing the entire sub-field)

    Log.debug("DetectorWithMappedPSF: resetting subfield array for new exposure.");

    reset();

    for (int binnumber=0; binnumber<wave_bins; binnumber++)  //%% Loop over different wavebins, for spectral dependency	
    {
            reset();

    // Integration (incl. jitter): point sources + background

    double centerRA, centerDec;  //%%
    tie(centerRA, centerDec) = camera.exposeDetector(*this, startTime, exposureTime, readoutTimeBeforeNextExposure, binnumber, subsubfieldx, subsubfieldy);	//%% added binnumber, subsubfield, ra/dec for spectral dependency

    // Convolve with the point spread function

    int fieldnumber = subsubfieldx * numsubsubfieldsy + subsubfieldy;  //%% Unique subfield ID to find correct vector position of PSF etc
    int psfnumber = binnumber + fieldnumber * wave_bins;  //%% unique psf number based on subfield and wavebin

    convolveWithPsf(psfnumber);  //%% added psfnumber for spectral dependency

    camera.SkyBackground(*this, startTime, exposureTime, readoutTimeBeforeNextExposure, binnumber, centerRA, centerDec, subsubfieldx, subsubfieldy);  //%% For spectral dependency - add the diffuse background separately after the convolution to eliminate edge effects

    // Apply flatfield (at sub-pixel level)

    if (includeFlatfield)
    {
        Log.debug("DetectorWithMappedPSF: applying Flatfield.");

        applyFlatfield(subsubfieldx, subsubfieldy);  //%% Added subsubfield for spectral dependence.
    }
    else
    {
        Log.debug("DetectorWithMappedPSF: no flatfield applied.");
    }

    // Rebin from a subpixel map to a pixel map

    Log.debug("DetectorWithMappedPSF: rebinning sub-pixel map into pixel map.");

    rebin();

    // Apply throughput efficiency on the pixel map
    // This takes into account the QE, vignetting, polarisation, and particulate & molecular contamination.
    // PixelMap units change from [photons] to [electrons]

    applyThroughputEfficiency(binnumber, subsubfieldx, subsubfieldy);

    // Brighter-Fatter effect

    if (includeBFE)
    {
        Log.debug("DetectorWithMappedPSF: adding Brighter-Fatter effect");

        applyBFE();
    }
    else
    {
        Log.debug("DetectorWithMappedPSF: no Brighter-Fatter effect added");
    }

    addSubSubField(subsubfieldx, subsubfieldy);  //%% Stitch the smaller pixel map additively to the arger complete map

}









/**
 * \brief: Add the given flux value to the value of the sub-pixel that corresponds to the given coordinates 
 *         in the focal plane. Return the pixel coordinates of the pixel to which the flux was added.
 *
 * \param xFP: X-coordinate of the sub-pixel in the focal plane in the FP reference frame [mm].
 * 
 * \param yFP: Y-coordinate of the sub-pixel in the focal plane in the FP reference frame [mm].
 * 
 * \param flux: Flux to add to the sub-pixel map [photons].
 *
 * \return (isInSubfield, row, col) 
 *         isInSubfield: True if (xFP, yFP) are on the subfield, false otherwise;
 *         row: sub-field (not CCD) row number of the pixel to which the flux was added;
 *         col: sub-field (not CCD) column number of the pixel to which the flux was added.
 */

tuple<bool, double, double> DetectorWithMappedPSF::addFlux(double xFP, double yFP, double flux, int subsubfieldx, int subsubfieldy) //%% Added subsubfield for spectral dependence
{
    // Convert from FP coordinates to CCD pixel coordinates

    double pixRow, pixColumn;
    tie(pixRow, pixColumn) = focalPlaneToPixelCoordinates(xFP, yFP);
    pixRow -= subFieldZeroPointRow + subsubfieldx * (numRowsPixelMap - 2 * overlapx);  //%% Take subsubfield into account
    pixColumn -= subFieldZeroPointColumn + subsubfieldy * (numColumnsPixelMap - 2 * overlapy);

    // Check if the star falls in the subfield. If not, don't add any flux, but simply return.

    if (!isInPixelMap(pixRow, pixColumn, subsubfieldx, subsubfieldy))  //%% Added subsubfield
    {
        return make_tuple(false, pixRow, pixColumn);
    }

    // Sub-field coordinates, taking into account the edge pixels
    // (subpixRow, subpixColumn) are the indices of the star in the subpixelMap. So they are not
    // sub-pixel coordinates in the CCD frame, but in the subfield reference frame.
    // (no longer rounded since the implementation of charge diffusion and jitter smoothing)

    const double subpixColumn = (pixColumn + numEdgePixels) * numSubPixelsPerPixel;
    const double subpixRow = (pixRow + numEdgePixels) * numSubPixelsPerPixel;

    // Add the flux to the sub-pixel map

    if (isInSubPixelMap(subpixRow, subpixColumn))
    {
        if (includeChargeDiffusion || includeJitterSmoothing)
        {
            // Apply either charge diffusion or jitter smoothing
            // - charge diffusion: kernel width -> configuration parameter "ChargeDiffusionStrength" [pixels]
            // - jitter smoothing: kernel width = 0.5 sub-pixels

            applyDiffusionKernel(subpixRow, subpixColumn, flux);
        }
        else
        {
            subPixelMap((int)floor(subpixRow), (int)floor(subpixColumn)) += flux;
        }

        return make_tuple(true, pixRow, pixColumn);
    }
    else
    {
        return make_tuple(false, pixRow, pixColumn);
    }
}

/**
 * \brief: Applies charge diffusion or jitter smoothing for the given flux at the given position in the sub-pixel map.
 *
 * \param subpixelRow: Row index [sub-pixels]. NOT a coordinate in the CCD frame, but in the subfield frame.
 * 
 * \param subpixColumn: Column index [sub-pixels].  NOT a coordinate in the CCD frame, but in the subfield frame.
 * 
 * \param flux: Flux for which to apply charge diffusion or jitter smoothing [photons].
 */
void DetectorWithMappedPSF::applyDiffusionKernel(double subpixRow, double subpixColumn, double flux)
{
    int sx = subpixColumn - (diffusionKernelImageSize - 1) / 2;
    int sy = subpixRow - (diffusionKernelImageSize - 1) / 2;

    double ox = subpixColumn - floor(subpixColumn);
    double oy = subpixRow - floor(subpixRow);

    // Establish diffusion kernel image

    signalResponse = IntegralOfAnalyticSignalResponse(diffusionKernelImageSize);
    signalResponse.addPart(ox, oy, 1., diffusionKernelWidth);

    for (unsigned int row = 0; row < diffusionKernelImageSize; row++)
    {
        for (unsigned int column = 0; column < diffusionKernelImageSize; column++)
        {
            diffusionKernel(row, column) = signalResponse(column, row); // Normalisation done by ()-operator
        }
    }

    // Add the flux to the sub-pixel map

    arma::span rowSpan = arma::span(max(0, sy), min((int)numRowsSubPixelMap, sy + diffusionKernelImageSize) - 1);
    arma::span columnSpan = arma::span(max(0, sx), min((int)numColumnsSubPixelMap, sx + diffusionKernelImageSize) - 1);

    arma::span xSpan = arma::span(max(0, sx) - sx, min((int)numColumnsSubPixelMap, sx + diffusionKernelImageSize) - sx - 1);
    arma::span ySpan = arma::span(max(0, sy) - sy, min((int)numRowsSubPixelMap, sy + diffusionKernelImageSize) - sy - 1);

    subPixelMap(rowSpan, columnSpan) += diffusionKernel(ySpan, xSpan) * flux;
}

/**
 * \brief Check whether the given (row, column) indices are within the array range of the subpixel map.
 *
 * \details The input parameters row & column come from a coordinate transformation
 *          in the focal plane, and as a result are not necessarily integers. For this 
 *          function it's not necessary to round them to the nearest integer. 
 *
 * \param row: Row index. NOT a coordinate in the CCD frame, but in the subfield frame. [sub-pixel].
 * 
 * \param column: Column index.NOT a coordinate in the CCD frame, but in the subfield frame.  [sub-pixel].
 *
 * \return  True if the given (row, column) coordinates are in the sub-pixel map; false otherwise.
 */
bool DetectorWithMappedPSF::isInSubPixelMap(double row, double column)
{
    return (column >= 0) && (row >= 0) && (column < numColumnsSubPixelMap) && (row < numRowsSubPixelMap);
}

/**
 * \brief: Add the given flux value to (all sub-pixels of) the sub-pixel map.
 *
 * \param flux: Flux to add to the sub-pixel map [photons/pixel].
 */

void DetectorWithMappedPSF::addFlux(double flux, int subsubfieldx, int subsubfieldy)  //%% For spectral dependency: only add flux to the core regions of the subfield, is extended if at the edge
{
    // The flux is expressed in [photons/pixel] but we need the quantity expressed
    // in [photons/subpixel]. There are (numSubPixelsPerPixel)^2 per pixel (the
    // name is thus a bit of a misnomer.).

    double fluxsub = flux / numSubPixelsPerPixel / numSubPixelsPerPixel;

    int colovlow = overlapy;
    int colovhigh = overlapy;
    int rowovlow = overlapx;
    int rowovhigh = overlapx;

    if (subsubfieldx == 0)
    {
        rowovlow = 0;
    }
    if (subsubfieldx == numsubsubfieldsx - 1)
    {
        rowovhigh = 0;
    }

    if (subsubfieldy == 0)
    {
        colovlow = 0;
    }
    if (subsubfieldy == numsubsubfieldsy - 1)
    {
        colovhigh = 0;
    }

    for (int i = numSubPixelsPerPixel * rowovlow; i < numRowsSubPixelMap - numSubPixelsPerPixel * rowovhigh; i++)
    {
        for (int j = numSubPixelsPerPixel * colovlow; j < numColumnsSubPixelMap - numSubPixelsPerPixel * colovhigh; j++)
        {
	    subPixelMap(i, j) += fluxsub;
        }
    }
}

/**
 * \brief: Multiply the sub-pixel map with the flatfield.
 * 
 * NOTE: The sub-pixel map contains extra edge pixels, but the flatfield
 *       map does not. These edge pixels are excluded from this flatfield
 *       multiplication.
 *
 * \pre Unit of the sub-pixels: [photons].
 * 
 * \pre Flatfield map at sub-pixel level, excl. edge pixels.
 * 
 * \pre Pixel, bias register, and smearing maps filled with zeroes.
 *
 * \post Pixel value in the sub-pixel map: [photons].
 * 
 * \post Pixel, bias, and smearing maps filled with zeroes.
 */

void DetectorWithMappedPSF::applyFlatfield(int subsubfieldx, int subsubfieldy)  //%% Added subsubfield for spectral dependence
{
    const unsigned int numEdgeSubPixels = numEdgePixels * numSubPixelsPerPixel;
    const unsigned int beginRow = numEdgeSubPixels;
    const unsigned int beginCol = numEdgeSubPixels;
    const unsigned int endRow = numRowsSubPixelMap - numEdgeSubPixels - 1;
    const unsigned int endCol = numColumnsSubPixelMap - numEdgeSubPixels - 1;

    const unsigned int FFbeginRow = numSubPixelsPerPixel * subsubfieldx * (numRowsPixelMap - 2 * overlapx) ;  //%% Only applying the FF of the corresponding subfield region
    const unsigned int FFbeginCol = numSubPixelsPerPixel * subsubfieldy * (numColumnsPixelMap - 2 * overlapy) ;
    const unsigned int FFendRow = FFbeginRow + numRowsSubPixelMap - numEdgeSubPixels - 1;
    const unsigned int FFendCol = FFbeginCol + numColumnsSubPixelMap - numEdgeSubPixels - 1;
    
    subPixelMap.submat(beginRow, beginCol, endRow, endCol) = subPixelMap.submat(beginRow, beginCol, endRow, endCol) % flatfieldMap.submat(FFbeginRow, FFbeginCol, FFendRow, FFendCol);
}

/**
 * \brief Rebin the sub-pixel map to pixel level and crop the edge pixels.
 *
 * \pre Unit of the pixel value in the sub-pixel map: [photons].
 * 
 * \pre Pixel, bias register, and smearing map filled with zeroes.
 *
 * \post Unit of pixel values in the sub-pixel map: [photons].
 * 
 * \post Bias register, and smearing maps filled with zeroes.
 */
void DetectorWithMappedPSF::rebin()
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
 * \brief      Set the Point Spread Function map for the subfield
 * 
 * \details    The PSF that is selected is dependent on the user input.
 */

void DetectorWithMappedPSF::setPsfForSubfield(int subsubfieldx, int subsubfieldy)
{
    // There is one PSF for the entire subfield, which we take the one of the center 
    // of the subfield.

    double xFPmm, yFPmm;

    tie(xFPmm, yFPmm) = getFocalPlaneCoordinatesOfSubfieldCenter(subsubfieldx, subsubfieldy);  //%% For spectral dependency, use correct subsubfield

    // Get the 'user specified' angular distance to the optical axis from the psf.
    // If the user didn't specify an angular distance, calculate it from the given
    // focal plane coordinates.

    double radius = psf->getRequestedDistanceToOpticalAxis();

    if (radius < 0.0)
    {
        radius = camera.getGnomonicRadialDistanceFromOpticalAxis(xFPmm, yFPmm);
    }

    int fieldnumber = subsubfieldx * numsubsubfieldsy + subsubfieldy; //%% unique number for subsubfield
    int fieldmax = numsubsubfieldsx * numsubsubfieldsy -1; //%% total number of subsubfields -1 for counting start at 0

    psf->select(radius, fieldnumber, fieldmax);  //%% added fieldnumber and max for spectral dependence


    // Get the 'user specified' orientation angle from the psf.
    // if the user didn't specify a rotation angle, calculate it 
    // from the given focal plane coordinates.

    double angle = psf->getRequestedRotationAngle();

    if (angle < 0.0)
    {
        angle = atan2(yFPmm, xFPmm);
    }

    //  Compensate for the orientation of the CCD wrt focal plane orientation.

    angle -= orientationAngle;
    psf->rotate(angle, fieldnumber, fieldmax);  //%% added fieldnumber and max for spectral dependence

    // Rebin the psfMap to the number of sub-pixels per pixel used for the Detector

    psfVector = psf->rebinToSubPixels(numSubPixelsPerPixel, fieldnumber); //%% fieldnumber to rebin the correct one

    // Allow the convolver to precompute some stuff given the PSF, so that it doesn't
    // need to be recomputed every convolution.

    for (int binnumber=0; binnumber<wave_bins; binnumber++)  //%% For spectral dependence: Loop over all wavebins and intialize all corresponding PSFs
    {
        convolver.initialise(numRowsSubPixelMap, numColumnsSubPixelMap, psfVector[binnumber], binnumber, wave_bins, fieldnumber, fieldmax);
    }

}












/**
 * \brief: Convolve the sub-pixel map with the PSF, keeping the same dimensions.
 *
 * \param psf: PSF.
 */

void DetectorWithMappedPSF::convolveWithPsf(int psfnumber)
{

    if (includeConvolution)
    {
        Log.debug("DetectorWithMappedPSF: convolving subPixelMap with PSF.");

        // subpixelMap serves here both as input as well as output matrix;

        convolver.convolve(subPixelMap, subPixelMap, psfnumber);  //%% Added psfnumber for spectral dependence
    }
    else
    {
        Log.debug("DetectorWithMappedPSF: no convolution applied.");
    }
}

/**
 * \brief: Creates the group(s) in the HDF5 file where the detector specific
 *         information will be stored.  These groups have to be created once,
 *         at the very beginning.
 */
void DetectorWithMappedPSF::initHDF5Groups()
{
    // Init the groups specific for the MappedPSF detector

    if (writeSubPixelImagesToHDF5)
    {
        hdf5File.createGroup("/SubPixelImages");
    }
}

/**
 * \brief: Writes the subpixel map for the HDF5 file.
 */

void DetectorWithMappedPSF::writeSubPixelMapToHDF5(int exposureNr)
{
    stringstream myStream;
    myStream << "subPixelImage" << setfill('0') << setw(6) << exposureNr;
    string imageName = myStream.str();

    // Add the image to the "SubPixelImages" group

    hdf5File.writeArray("/SubPixelImages", imageName, subPixelMap);
}
