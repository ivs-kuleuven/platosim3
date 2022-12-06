#include "DetectorWithMappedPSF.h"



/**
 * \brief Constructor.
 *
 * \details
 *
 * The constructor initializes the groups in the HDF5 file where the different maps (i.e. pixel map,
 * bias register map, smearing map, etc.) will be saved.
 *
 * The following maps are initialized to zero (partly through the base class Detector):
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
 * \param readoutTimeBeforeNextExposure Duration of the readout that takes place before the next exposure can start.
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
    flatfieldMap.ones(numRowsSubPixelMap, numColumnsSubPixelMap);

    if (includeFlatfield)
    {
        // Generate the flatfield map

        generateFlatfieldMap();
    }

    // Initialize and load the PSF. This will open the PSF HDF5 file and perform some basic checking,
    // Then select the proper PSF for the given subfield. Should only be done after calling configure().

    psf = new PointSpreadFunction(configParam, hdf5file);
    setPsfForSubfield();
}




/**
 * Destructor.
 */
DetectorWithMappedPSF::~DetectorWithMappedPSF()
{
    flushOutput();
    delete psf;
}


/**
 * \brief Configure the DetectorWithMappedPSF object using the given
 *        configuration parameters.
 *
 * \param configParam: Configuration parameters.
 */
void DetectorWithMappedPSF::configure(ConfigurationParameters &configParam)
{
    flatfieldNoiseRMS = configParam.getDouble("CCD/FlatfieldNoiseRMS");
    includeFlatfield = configParam.getBoolean("CCD/IncludeFlatfield");
    includeConvolution = configParam.getBoolean("CCD/IncludeConvolution");

    writeSubPixelImagesToHDF5 = configParam.getBoolean("ControlHDF5Content/WriteSubPixelImages");

    numSubPixelsPerPixel = configParam.getInteger("SubField/SubPixels");

    // Configuration parameters for the noise source random seeds

    flatfieldSeed = configParam.getLong("RandomSeeds/FlatFieldSeed");

    // Treat the specific configurations for a Mapped PSF

    string psfModel = configParam.getString("PSF/Model");

    if (psfModel == "MappedFromFile")
    {
        includeChargeDiffusion = configParam.getBoolean("PSF/" + psfModel + "/IncludeChargeDiffusion");
        includeJitterSmoothing = configParam.getBoolean("PSF/" + psfModel + "/IncludeJitterSmoothing");
        chargeDiffusionStrength = configParam.getDouble("PSF/" + psfModel + "/ChargeDiffusionStrength");

        if (includeChargeDiffusion)
        {
            generateDiffusionKernel(chargeDiffusionStrength * numSubPixelsPerPixel);
        }
        else if (includeJitterSmoothing)
        {
            generateDiffusionKernel(0.5);
        }
    }

    // Derive the dimensions of the sub-pixel map

    numRowsSubPixelMap = numRowsPixelMap * numSubPixelsPerPixel;       // TODO Add edge pixels
    numColumnsSubPixelMap = numColumnsPixelMap * numSubPixelsPerPixel; // TODO Add edge pixels

    // The configuration for the HDF5 contents

    writeFlatfieldMap = configParam.getBoolean("ControlHDF5Content/WriteFlatfieldMap");
    writeDiffusedPSF = configParam.getBoolean("ControlHDF5Content/WriteDiffusedPSF");
}






/**
 * \brief Set the PSF map for the sub-field and sets  the distortion map
 *
 * \details The PSF that is selected is dependent on the user input.
 */
void DetectorWithMappedPSF::setPsfForSubfield()
{
    // There is one PSF for the entire subfield, which we take the one of the center
    // of the subfield.

    double xFPmm, yFPmm;
    tie(xFPmm, yFPmm) = getFocalPlaneCoordinatesOfSubfieldCenter();

    psf->select(xFPmm, yFPmm);
    distortionMap = psf->getDistortionMap();

    if(psf->getNumSubPixelsPerPixel() < numSubPixelsPerPixel)
    {
        throw IllegalArgumentException(string("DetectorWithMappedPSF.setPsfForSubfield: ") + 
            "The sub-pixel resolution of the PSF (" + to_string(psf->getNumSubPixelsPerPixel()) +
                    ") must be at least that of the sub-field (" + to_string(numSubPixelsPerPixel) + ")");
    }

    // If requestied save the diffused PSF to the output file
    if (writeDiffusedPSF)
    {
      writeDiffusedPSFToHDF5(psf);
    }

    //  Compensate for the orientation of the CCD wrt focal plane orientation.

    psf->rotate(-rotationAnglePsf);

    // Rebin the psfMap to the number of sub-pixels per pixel used for the Detector

    psfMap = psf->rebinToSubPixels(numSubPixelsPerPixel);

    // Allow the convolver to precompute some stuff given the PSF, so that it doesn't
    // need to be recomputed every convolution.

    convolver.initialise(numRowsSubPixelMap, numColumnsSubPixelMap, psfMap);
}







/**
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

    unsigned int Nrows = 2 * numRowsPixelMap * numSubPixelsPerPixel;
    unsigned int Ncolumns = 2 * numColumnsPixelMap * numSubPixelsPerPixel;

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

    arma::Mat<float> prnu(numRowsPixelMap, numColumnsPixelMap, arma::fill::zeros);

    for (unsigned int row = 0; row < numRowsPixelMap; row++)
    {
        for (unsigned int column = 0; column < numColumnsPixelMap; column++)
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
 *         This differs from the normal Detector::reset() because it includes resetting the subPixelMap.
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
    // Advance the internal clock until the given start time

    internalTime = startTime;

    // Clear all arrays

    Log.debug("Detector: resetting subfield array for new exposure.");
    reset();

    // Integration of point sources and background, taking into account jitter + drift.

    Log.info("DetectorWithMappedPSF: Integrating light for exposure " + to_string(exposureNr) + " with exposure time = " + to_string(exposureTime));

    integrateLight(exposureNr, startTime, exposureTime);

    // If this is the first exposure, we should initialize the number of occupied traps.
    // This can only be done after the detector
    // has been exposed to the skybackground.
    // => Check if CTI is included && We use the Short2013 model

    if (exposureNr == beginExposureNr) {
      if (includeCTIeffects &&
          (CTImodel == "Short2013" || CTImodel == "Short2013FromFile"))
      {
          setInitialNumberOfOccupiedTraps(numberOfOccupiedTrapsPixelMap);
      }
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
    // Write the cosmic hits to the HDF5 file

        Log.debug("Detector: Writing Cosmics of the PixelMap, smearing map, bias map #" + to_string(exposureNr) + " to HDF5 file.");

    writeCosmicHitsToHDF5(exposureNr);

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
void DetectorWithMappedPSF::integrateLight(int exposureNr, double startTime, double exposureTime)
{
    // Reset the sub-field (i.e. get rid of the previous exposure, by zeroing the entire sub-field)

    Log.debug("DetectorWithMappedPSF: resetting subfield array for new exposure.");

    reset();

    // Integration (incl. jitter): point sources

    camera.exposeDetectorWithStars(*this, startTime, exposureTime, readoutTimeBeforeNextExposure);

    // Convolve with the point spread function

    convolveWithPsf();

    // Integration: background

    camera.exposeDetectorWithSkyBackground(*this, startTime, exposureTime, readoutTimeBeforeNextExposure);

    // Apply flatfield (at sub-pixel level)

    if (includeFlatfield)
    {
        Log.debug("DetectorWithMappedPSF: applying Flatfield.");

        applyFlatfield();
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

    applyThroughputEfficiency();

    // Apply the charge injection which will mitigate the CTI. The injection happens in electrons, 
    // so the throughput efficiency should already have been applied. In principle, the injected charges do 
    // feel the PRNU, but for the MappedPSF we first need to apply the PRNU on sub-pixel level and afterwards
    // apply the throughputEfficiency() at pixel level, so there is no possibilty to respect the order
    // (1) throughput (2) charge injection (3) PRNU.

    if (includeChargeInjection)
    {
        Log.debug("Detector: applying charge injection");
        applyChargeInjection();
    }



    // Apply the effects of readout smearing due to an open shutter. Because there is no shutter,
    // the pixels are still receiving photons from the sky, while they are being transfered towards
    // the readout register.
    // Pixel units before: [electrons]
    // Pixel units after: [electrons]


    if (includeOpenShutterSmearing)
    {
        Log.debug("Detector: applying open shutter smearing.");
        applyOpenShutterSmearing(exposureTime);
    }
    else
    {
         Log.debug("Detector: no open shutter smearing applied.");
    }

    // Apply poisson distributed photon noise
    // Pixel units before: [electrons]
    // Pixel units after: [electrons]


    if (includePhotonNoise)
    {
        Log.debug("Detector: adding photon noise.");
        addPhotonNoise();
    }
    else
    {
        Log.debug("Detector: no photon noise added.");
    }


    // Add dark current

    if (includeDarkSignal)
    {
        Log.debug("DetectorWithMappedPSF: adding dark current");

        addDarkSignal(exposureTime);
    }
    else
    {
        Log.debug("Detector: no dark current added");
    }

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

tuple<bool, double, double> DetectorWithMappedPSF::addFlux(double xFP, double yFP, double flux)
{
    // Convert from FP coordinates to real-valued CCD pixel coordinates

    double pixRow, pixColumn;
    tie(pixRow, pixColumn) = focalPlaneToPixelCoordinates(xFP, yFP);
    pixRow -= subFieldZeroPointRow;
    pixColumn -= subFieldZeroPointColumn;

    // Check if the star falls in the subfield. If not, don't add any flux, but simply return.

    if (!isInPixelMap(pixRow, pixColumn))
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
 * \brief Insert the extended ghost with the given radius and flux at the given focal-plane position.
 * 
 * Note that the extended source will be convolved with the PSF in a next step.
 *
 * \param x0: Focal-plane x-coordinate of the centre of the extended ghost [mm].
 * \param y0: Focal-plane y-coordinate of the centre of the extended ghost [mm].
 * \param radius: Radius of the extended ghost [mm].
 * \param flux: Flux of the extended ghost [photons].
 *
 * \return: Whether or not the extended source falls (at least partially) on the sub-field, and the
 *          (row, column) coordinates of the centre of the extended ghost in the pixel map.
 */
tuple<bool, double, double> DetectorWithMappedPSF::addExtendedGhost(double x0, double y0, double radius, double flux)
{
    // Calculate the number of sub-pixels in the extended ghost

    double radiusSubPixels = radius * 1000 / pixelSize * numSubPixelsPerPixel;    // Radius [sub-pixels]
    double radiusSubPixelsSquared = pow(radiusSubPixels, 2);                      // Squared radius [sub-pixels^2]

    double numSubPixels = PI * pow(radiusSubPixels, 2);       // Area of the extended ghost [sub-pixels]
    double fluxPerSubPixel = flux / numSubPixels;             // Flux [photons / sub-pixel]

    // Calculate the (row, column) coordinates of the centre of the extended source in the sub-pixel map

    double row0, column0;
    tie(row0, column0) = focalPlaneToPixelCoordinates(x0, y0);
    row0 -= subFieldZeroPointRow;
    column0 -= subFieldZeroPointColumn;

    row0 *= numSubPixelsPerPixel;
    column0 *= numSubPixelsPerPixel;

    bool ghostInSubPixelMap = false;

    // Try to add flux to all pixels covered by the extended ghosts

    for(int row = row0 - radiusSubPixels; row <= row0 + radiusSubPixels; row++)
    {
        for(int column = column0 - radiusSubPixels; column <= column0 + radiusSubPixels; column++)
        {
            if (isInSubPixelMap(row, column) && pow(column - column0, 2) + pow(row - row0, 2) <= radiusSubPixelsSquared)
            {
                ghostInSubPixelMap = true;
                subPixelMap(row, column) += fluxPerSubPixel;
            }
        }
    }

    return  make_tuple(ghostInSubPixelMap, row0 / numSubPixelsPerPixel, column0 / numSubPixelsPerPixel);
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
    int coveredSubPixelsLeft   = coveredLeft * numSubPixelsPerPixel;
    int coveredSubPixelsRight  = coveredRight * numSubPixelsPerPixel;
    int coveredSubPixelsBottom = coveredBottom * numSubPixelsPerPixel;
    int coveredSubPixelsTop     = coveredTop * numSubPixelsPerPixel;

    return (column >= coveredSubPixelsLeft) && (row >= coveredSubPixelsBottom) && (column < numColumnsPixelMap*numSubPixelsPerPixel - coveredSubPixelsRight) && (row < numRowsPixelMap*numSubPixelsPerPixel - coveredSubPixelsTop);

}

/**
 * \brief: Add the given flux value to (all sub-pixels that are not covered by a metallic
 *         shield of) the sub-pixel map.
 *
 * \param flux: Flux to add to the sub-pixel map [photons/pixel].
 */
void DetectorWithMappedPSF::addFlux(double flux)
{
    int coveredSubPixelsLeft   = coveredLeft * numSubPixelsPerPixel;
    int coveredSubPixelsRight  = coveredRight * numSubPixelsPerPixel;
    int coveredSubPixelsBottom = coveredBottom * numSubPixelsPerPixel;
    int coveredSubPixelsTop    = coveredTop * numSubPixelsPerPixel;

    bool isBlockedOff = (coveredSubPixelsLeft + coveredSubPixelsRight >= numColumnsSubPixelMap || coveredSubPixelsBottom + coveredSubPixelsTop >= numRowsSubPixelMap);

    if (!isBlockedOff)
    {
      subPixelMap.submat(coveredSubPixelsBottom, coveredSubPixelsLeft,
                         numRowsSubPixelMap - coveredSubPixelsTop - 1,
                         numColumnsSubPixelMap - coveredSubPixelsRight - 1) +=
          flux / numSubPixelsPerPixel / numSubPixelsPerPixel;
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
void DetectorWithMappedPSF::applyFlatfield()
{
    const unsigned int numEdgeSubPixels = numEdgePixels * numSubPixelsPerPixel;
    const unsigned int beginRow = numEdgeSubPixels;
    const unsigned int beginCol = numEdgeSubPixels;
    const unsigned int endRow = numRowsSubPixelMap - numEdgeSubPixels - 1;
    const unsigned int endCol = numColumnsSubPixelMap - numEdgeSubPixels - 1;

    subPixelMap.submat(beginRow, beginCol, endRow, endCol) = subPixelMap.submat(beginRow, beginCol, endRow, endCol) % flatfieldMap;
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
 * \brief: Convolve the sub-pixel map with the PSF, keeping the same dimensions.
 *
 * \param psf: PSF.
 */
void DetectorWithMappedPSF::convolveWithPsf()
{

    if (includeConvolution)
    {
        Log.debug("DetectorWithMappedPSF: convolving subPixelMap with PSF.");

        // subpixelMap serves here both as input as well as output matrix;

        convolver.convolve(subPixelMap, subPixelMap);
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







/**
 * \brief: Write the diffused and rotated PSF to the HDF5 file.
 */

void DetectorWithMappedPSF::writeDiffusedPSFToHDF5(PointSpreadFunction *psf)
{
    arma::fmat psfMap = psf->getOriginalPSF();
    arma::fmat diffusedPsf = arma::fmat(size(psfMap), arma::fill::zeros);

    int numRows = size(psfMap)(0);
    int numColumns = size(psfMap)(1);

    // set the diffusion kernel image size
    int psfSubPixelsPerPixel = psf->getNumSubPixelsPerPixel();
    generateDiffusionKernel(chargeDiffusionStrength*psfSubPixelsPerPixel);

    for (int row=0; row < numRows; row++)
    {
      for (int column=0; column < numColumns; column++)
      {
                applyDiffusionKernelOnPSF(row, column, psfMap(row, column), diffusedPsf, psfSubPixelsPerPixel);
      }
    }

    // reset the diffusion kernel
    generateDiffusionKernel(chargeDiffusionStrength*numSubPixelsPerPixel);


    // rotate the diffused PSF
    diffusedPsf = ArrayOperations::rotateArray(diffusedPsf, -rotationAnglePsf);
    diffusedPsf /= arma::accu(diffusedPsf);

    // write the diffused psf to the output hdf5 file
    hdf5File.writeArray("/PSF", "diffusedPSF", diffusedPsf);
}







void DetectorWithMappedPSF::applyDiffusionKernelOnPSF(double subpixRow, double subpixColumn, double flux, arma::fmat& psf, int numberOfPsfSubpixelsPerPixel)
{
    int sx = subpixColumn - (diffusionKernelImageSize - 1) / 2;
    int sy = subpixRow - (diffusionKernelImageSize - 1) / 2;

    double ox = subpixColumn - floor(subpixColumn);
    double oy = subpixRow - floor(subpixRow);

    int numRows = size(psf)(0);
    int numColumns = size(psf)(1);

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

    // Add the flux to the psf


    arma::span rowSpan = arma::span(max(0, sy), min((int)numRows, sy + diffusionKernelImageSize) - 1);
    arma::span columnSpan = arma::span(max(0, sx), min((int)numColumns, sx + diffusionKernelImageSize) - 1);

    arma::span xSpan = arma::span(max(0, sx) - sx, min((int)numColumns, sx + diffusionKernelImageSize) - sx - 1);
    arma::span ySpan = arma::span(max(0, sy) - sy, min((int)numRows, sy + diffusionKernelImageSize) - sy - 1);

    psf(rowSpan, columnSpan) += diffusionKernel(ySpan, xSpan) * flux;
}


/*
 * /brief: determins if three points in (given in an array) are colinear
 * /input: an array with the points of the form {{x1, y1}, {x2, y2}, {x3, y3}}
 * /note:
 * Three points are colinear if the determinant of the matrix
 * | 1  1  1|
 * |x1 x2 x3| is equal to zero.
 * |y1 y2 y3|
 */
bool DetectorWithMappedPSF::areColinear(std::array<std::array<double, 2>, 3> points)
{
    double determinant = points[1][0] * points[2][1] - points[2][0] * points[1][1] - points[0][0] * points[2][1] + points[2][0] * points[0][1] + points[0][0] * points[1][1] - points[1][0] * points[0][1];
    if(abs(determinant) < 0.01)
    {
      return true;
    }
    else
    {
      return false;
    }
}




/*
 * /brief: applies the field distortion on the inputparameters from the distortion map
 * /input: FP coordinates [mm]
 */
void DetectorWithMappedPSF::applyDistortion(double &x, double &y)
{
    // We try to select 4 undistorted "closest"* points around the input point
    // together with their respectuve distorted points. To know
    // wether a point should be included instead of another we compare their distance*
    // the 4 points with the smallest distance* get selected.

    // * REMARK: The distance we use is not with the euclidean metric. This is because we
    // want to select 4 points that describe a square in which the input point lies.

    std::array<std::array<double, 2>,4> ClosestUndistortedCoordinates;
    std::array<std::array<double, 2>,4> ClosestDistortedCoordinates;

    std::array<double,4> minDistance;
    minDistance.fill(std::numeric_limits<double>::max());

    for (auto& coordinates : distortionMap)
    {
        double distance = std::max( pow(coordinates[0] - x, 2), pow(coordinates[1] - y, 2));

        // for a point to be included, its distance to the input point should be at least smaller
        // then the highest distance that we keep track of.

        if(distance < minDistance[3])
        {
            // Once we know we want to include the point, we check where its position in the list should
            // be so that the order nearest -> furthest remains.
            int j=2;
            while( (j >= 0) && (distance < minDistance[j]))
            {
                minDistance[j+1] = minDistance[j];
                ClosestUndistortedCoordinates[j+1] = ClosestUndistortedCoordinates[j];
                ClosestDistortedCoordinates[j+1] = ClosestDistortedCoordinates[j];
                j = j-1;
            }
            minDistance[j+1] = distance;
            ClosestUndistortedCoordinates[j+1] = {coordinates[0], coordinates[1]};
            ClosestDistortedCoordinates[j+1] = {coordinates[2], coordinates[3]};
        }
    }

    // We now have selected 4 points that form a square in which the input point will lie. We now want to
    // order these points so that the lowest, most left point corresponds to index 0, and the next indices
    // follow the points counterclockwise.

    for (int i=0; i<3; i++)
    {
        bool isBigger   = (ClosestUndistortedCoordinates[i])[0] > (ClosestUndistortedCoordinates[i+1])[0];

        if(isBigger)
        {
            for (int k=i; k>=0; k--)
            {
                std::array<double, 2> dummy        = ClosestUndistortedCoordinates[k];
                ClosestUndistortedCoordinates[k]   = ClosestUndistortedCoordinates[k+1];
                ClosestUndistortedCoordinates[k+1] = dummy;
                dummy                              = ClosestDistortedCoordinates[k];
                ClosestDistortedCoordinates[k]     = ClosestDistortedCoordinates[k+1];
                ClosestDistortedCoordinates[k+1]   = dummy;
            }

        }
    }

  for (int i=0; i<3; i++)
  {
      bool isBigger = (ClosestUndistortedCoordinates[i])[0] < (ClosestUndistortedCoordinates[i+1])[0];
      bool secondIsBigger = (ClosestUndistortedCoordinates[i])[1] > (ClosestUndistortedCoordinates[i+1])[1];

      if(!isBigger && secondIsBigger)
      {
          std::array<double, 2> dummy        = ClosestUndistortedCoordinates[i];
          ClosestUndistortedCoordinates[i]   = ClosestUndistortedCoordinates[i+1];
          ClosestUndistortedCoordinates[i+1] = dummy;
          dummy                              = ClosestDistortedCoordinates[i];
          ClosestDistortedCoordinates[i]     = ClosestDistortedCoordinates[i+1];
          ClosestDistortedCoordinates[i+1]   = dummy;
      }
  }

  // We approximate the distorted coordinates using a linear combination of the 4 distorted coordinates corresponding
  // with the points around the input point.

  std::array<double, 4> constants;
  std::array<int, 4> oppositeIdx = {3, 2, 1, 0};

  // We get the constants of the linear combination
  double area = (ClosestUndistortedCoordinates[0][0] - ClosestUndistortedCoordinates[3][0]) * (ClosestUndistortedCoordinates[0][1] - ClosestUndistortedCoordinates[3][1]);
  for (int i=0; i<4; i++)
  {
      std::array<double, 2> oppositePoint = ClosestUndistortedCoordinates[oppositeIdx[i]];
      constants[i] = abs( (oppositePoint[0]-x)*(oppositePoint[1]-y) ) / area;
  }

  x = constants[0] * ClosestDistortedCoordinates[0][0] + constants[1] * ClosestDistortedCoordinates[1][0] + constants[2] * ClosestDistortedCoordinates[2][0] + constants[3] * ClosestDistortedCoordinates[3][0];
  y = constants[0] * ClosestDistortedCoordinates[0][1] + constants[1] * ClosestDistortedCoordinates[1][1] + constants[2] * ClosestDistortedCoordinates[2][1] + constants[3] * ClosestDistortedCoordinates[3][1];

  // Round up/down the numbers to 4 decimal numbers to get consistent results

  float x_dummy = (int) (x * 10000 + 0.5);
  float y_dummy = (int) (y * 10000 + 0.5);
  x = (float) x_dummy / 10000;
  y = (float) y_dummy / 10000;
}


/*
 * /brief: applies the inverse of the field distortion on the input coordinates
 * /input: FP coordinates [mm]
 */
void DetectorWithMappedPSF::applyInverseDistortion(double &x, double &y)
{
    // We try to select 4 distorted "closest"* points around the input point
    // together with their respectuve undistorted points. To know
    // wether a point should be included instead of another we compare their distance*
    // the 4 points with the smallest distance* get selected.

    // * REMARK: The distance we use is not with the euclidean metric. This is because we
    // want to select 4 points that describe a square in which the input point lies.

    std::array<std::array<double, 2>,4> ClosestUndistortedCoordinates;
    std::array<std::array<double, 2>,4> ClosestDistortedCoordinates;

    std::array<double,4> minDistance;
    minDistance.fill(std::numeric_limits<double>::max());

    for (auto& coordinates : distortionMap)
    {
        double distance = std::max( pow(coordinates[2] - x, 2), pow(coordinates[3] - y, 2));

        // for a point to be included, its distance to the input point should be at least smaller
        // then the highest distance that we keep track of.

        if(distance < minDistance[3])
        {
            // Once we know we want to include the point, we check where its position in the list should
            // be so that the order nearest -> furthest remains.
            int j=2;
            while( (j >= 0) && (distance < minDistance[j]))
            {
                minDistance[j+1] = minDistance[j];
                ClosestUndistortedCoordinates[j+1] = ClosestUndistortedCoordinates[j];
                ClosestDistortedCoordinates[j+1] = ClosestDistortedCoordinates[j];
                j = j-1;
            }
            minDistance[j+1] = distance;
            ClosestUndistortedCoordinates[j+1] = {coordinates[0], coordinates[1]};
            ClosestDistortedCoordinates[j+1] = {coordinates[2], coordinates[3]};
        }
    }

    // We now have selected 4 points that form a square in which the input point will lie. We now want to
    // order these points so that the lowest, most left point corresponds to index 0, and the next indices
    // follow the points counterclockwise.

    for (int i=0; i<3; i++)
    {
        bool isBigger   = (ClosestDistortedCoordinates[i])[0] > (ClosestDistortedCoordinates[i+1])[0];

        if(isBigger)
        {
            for (int k=i; k>=0; k--)
            {
                std::array<double, 2> dummy        = ClosestUndistortedCoordinates[k];
                ClosestUndistortedCoordinates[k]   = ClosestUndistortedCoordinates[k+1];
                ClosestUndistortedCoordinates[k+1] = dummy;
                dummy                              = ClosestDistortedCoordinates[k];
                ClosestDistortedCoordinates[k]     = ClosestDistortedCoordinates[k+1];
                ClosestDistortedCoordinates[k+1]   = dummy;
            }
        }
    }

  for (int i=0; i<3; i++)
  {
      bool isBigger = (ClosestDistortedCoordinates[i][0] < ClosestDistortedCoordinates[i+1][0]) && (abs(ClosestDistortedCoordinates[i][0] - ClosestDistortedCoordinates[i+1][0]) > 0.01) ;
      bool secondIsBigger = (ClosestDistortedCoordinates[i])[1] > (ClosestDistortedCoordinates[i+1])[1];

      if(!isBigger && secondIsBigger)
      {
          std::array<double, 2> dummy        = ClosestUndistortedCoordinates[i];
          ClosestUndistortedCoordinates[i]   = ClosestUndistortedCoordinates[i+1];
          ClosestUndistortedCoordinates[i+1] = dummy;
          dummy                              = ClosestDistortedCoordinates[i];
          ClosestDistortedCoordinates[i]     = ClosestDistortedCoordinates[i+1];
          ClosestDistortedCoordinates[i+1]   = dummy;
      }
  }

  // The distorted coordinates have been approximated as a linear combination of the
  // 4 nearest points. We can use this to know the coordinates of the undistorted point.

  std::array<double, 4> constants;
  std::array<int, 4> oppositeIdx = {3, 2, 1, 0};

  float sum = 0;
  for (int i=0; i<4; i++)
  {
      constants[i] = abs( (ClosestDistortedCoordinates[oppositeIdx[i]][0] - x)*(ClosestDistortedCoordinates[oppositeIdx[i]][1]-y));
      sum = sum + abs( (ClosestDistortedCoordinates[oppositeIdx[i]][0] - x)*(ClosestDistortedCoordinates[oppositeIdx[i]][1]-y));
  }
  for (int i=0; i<4; i++)
  {
      constants[i] = constants[i]/sum;
  }

  x = constants[0] * ClosestUndistortedCoordinates[0][0] + constants[1] * ClosestUndistortedCoordinates[1][0] + constants[2] * ClosestUndistortedCoordinates[2][0] + constants[3] * ClosestUndistortedCoordinates[3][0];
  y = constants[0] * ClosestUndistortedCoordinates[0][1] + constants[1] * ClosestUndistortedCoordinates[1][1] + constants[2] * ClosestUndistortedCoordinates[2][1] + constants[3] * ClosestUndistortedCoordinates[3][1];

  // Round up/down the numbers to 4 decimal numbers to get consistent results

  float x_dummy = (int) (x * 10000 + 0.5);
  float y_dummy = (int) (y * 10000 + 0.5);
  x = (float) x_dummy / 10000;
  y = (float) y_dummy / 10000;


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

void DetectorWithMappedPSF::generateThroughputMap()
{
    Log.info("DetectorWithMappedPSF: generating throughput map.");

    throughputMap.fill(1.0);

    if(includeRelativeTransmissivity  && includeOpenShutterSmearing)
        mechanicalVignettingMask.fill(1);

    double xFPmmDistorted, yFPmmDistorted;             // Distorted focal plan coordinates   [mm]
    double xFPmmUndistorted, yFPmmUndistorted;         // Undistorted focal plan coordinates [mm]
    double angle;                                      // Gnomonic radial distance from the optical axis [rad]
    double relativeTransmissivityVariation;


//    const double refAnglePolarizationRadians = deg2rad(refAnglePolarization);       // Reference angle for the polarisation efficiency [radians]
//    const double acosPolarizationEfficiency = acos(polarizationEfficiency);

//    const double refAngleQuantumEfficiencyRadians = deg2rad(refAngleQE);     // Reference angle for the quantum efficiency [radians]
//    const double acosQuantumEfficiency = acos(relativeRefEfficiencyQE);        // Relative efficiency due to the angle dependency of the QE at the reference angle

    if (includeRelativeTransmissivity || includePolarization || includeQuantumEfficiency)
    {
        // Loop over all pixels in the pixel map

        for (unsigned int row = 0; row < numRowsPixelMap; row++)
        {
            for (unsigned int column = 0; column < numColumnsPixelMap; column++)
            {
                // Distorted pixel coordinates (in the detector) -> distorted focal-plane coordinates

                tie(xFPmmDistorted, yFPmmDistorted) = pixelToFocalPlaneCoordinates(row + subFieldZeroPointRow, column + subFieldZeroPointColumn);
                xFPmmUndistorted = xFPmmDistorted;
                yFPmmUndistorted = yFPmmDistorted;
                // Convert from distorted to undistorted focal plane coordinates (Cf GitHub issue #716)

                applyInverseDistortion(xFPmmUndistorted, yFPmmUndistorted);

                // Angular distance [radians] of the pixel from the optical axis

                angle = camera.getGnomonicRadialDistanceFromOpticalAxis(xFPmmUndistorted, yFPmmUndistorted);  // [radians]

                if (includeRelativeTransmissivity)
                {
                    if (angle >= radiusFOV)
                    {
                        throughputMap(row, column) = 0.0;

                        if (includeOpenShutterSmearing)
                            mechanicalVignettingMask(row, column) = 0;
                    }

                    else
                    {
                        angle = rad2deg(angle); // [degrees]
                        relativeTransmissivityVariation = (  relTransmissivityCoefVector[0] * pow(angle, 2)
                                                           + relTransmissivityCoefVector[1] * pow(angle, 4)
                                                           + relTransmissivityCoefVector[2] * pow(angle, 6)) / 100.;

                        throughputMap(row, column) *= (1 - relativeTransmissivityVariation);
                    }
                }

                // Polarisation (Eq. 4-11 in PLATO-DLR-PL-RP-001)

                // NOTE: the polarization is angle dependent, but since no info on this dependency is currently available,
                //       we assume fow now it is fixed over the entire FOV.

                if (includePolarization)
                    throughputMap(row, column) *= expectedValuePolarization; //cos(angle / refAnglePolarizationRadians * acosPolarizationEfficiency);

                // Quantum efficiency (Eq. 4-12 in PLATO-DLR-PL-RP-001)
                // Pixel units before: [photons]
                // Pixel units after: [electrons]

                // NOTE: the QE is angle dependent, but since no info on this dependency is currently available,
                //       we assume for now it is fixed over the entire FOV.

                if (includeQuantumEfficiency)
                    throughputMap(row, column) *= meanQE * meanAngleDependencyQE; //(meanQE * cos(angle / refAngleQuantumEfficiencyRadians * acosQuantumEfficiency));
            }
        }
    }

    // Particulate contamination (Sect. 4.2.4.3 in PLATO-DLR-PL-RP-001)

    if (includeParticulateContamination)
    {
        throughputMap *= particulateContaminationEfficiency;
    }

    // Molecular contamination (Sect. 4.2.4.4 in PLATO-DLR-PL-RP-001)

    if (includeMolecularContamination)
    {
        throughputMap *= molecularContaminationEfficiency;
    }
}
