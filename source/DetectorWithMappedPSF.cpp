#include "DetectorWithMappedPSF.h"

/**
 * \brief Default constructor.
 * 
 * \param configParam: Configuration parameters.
 * 
 * \param hdf5file: HDF5 file from which to read the PSFs.
 */
DetectorWithMappedPSF::DetectorWithMappedPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure) : Detector(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure),
                                                                                                                                                                                                                                                                                                        includeFlatfield(true),
                                                                                                                                                                                                                                                                                                        writeSubPixelImagesToHDF5(false) {}

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

    int Nrows = 2 * numRowsPixelMap * numSubPixelsPerPixel;
    int Ncolumns = 2 * numColumnsPixelMap * numSubPixelsPerPixel;

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
    // Convert from FP coordinates to CCD pixel coordinates

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
    return (column >= 0) && (row >= 0) && (column < numColumnsSubPixelMap) && (row < numRowsSubPixelMap);
}

/**
 * \brief: Add the given flux value to (all sub-pixels of) the sub-pixel map.
 *
 * \param flux: Flux to add to the sub-pixel map [photons/pixel].
 */
void DetectorWithMappedPSF::addFlux(double flux)
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
