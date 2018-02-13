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
 */

DetectorWithMappedPSF::DetectorWithMappedPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator)
: Detector(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator), includeFlatfield(true), writeSubPixelImagesToHDF5(false)
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

    if(includeFlatfield)
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
    // Treat the specific configurations for a Mapped PSF

    flatfieldNoiseAmplitude   = configParam.getDouble("CCD/FlatfieldPtPNoise");
    includeFlatfield          = configParam.getBoolean("CCD/IncludeFlatfield");
    includeConvolution        = configParam.getBoolean("CCD/IncludeConvolution");

    writeSubPixelImagesToHDF5 = configParam.getBoolean("ControlHDF5Content/WriteSubPixelImages");

    numSubPixelsPerPixel    = configParam.getInteger("SubField/SubPixels");

    // Configuration parameters for the noise source random seeds

    flatfieldSeed           = configParam.getLong("RandomSeeds/FlatFieldSeed");

    // Derive the dimensions of the sub-pixel map

    numRowsSubPixelMap    = numRowsPixelMap    * numSubPixelsPerPixel;  // TODO Add edge pixels
    numColumnsSubPixelMap = numColumnsPixelMap * numSubPixelsPerPixel;  // TODO Add edge pixels

 }










 /**
 * \brief: Generate the (random) flatfield variations.  This map is generated
 *         at sub-pixel level but without the edge pixels.
 *
 * https://github.com/python-acoustics/python-acoustics/blob/master/acoustics/generator.py#L108
 */

void DetectorWithMappedPSF::generateFlatfieldMap()
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
    flatfieldMap *= flatfieldNoiseAmplitude;    // [0, flatfialdNoiseAmplitude]
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
 * \brief: Zeroes the pixel, bias register, and the smearing maps.
 *
 * \pre pixel, bias register, and smearing maps filled with values from previous exposure.
 *
 * \post pixel, bias register, and smearing maps filled with zeroes.
 */

void DetectorWithMappedPSF::reset()
{
    pixelMap.zeros();
    biasMap.zeros();
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
 * \param exposureNr:   Sequential number of the exposure
 * \param startTime:    Starting time of the exposure [s].
 * \param exposureTime: Duration of the exposure [s].
 * 
 * \return endTime:     Time after the exposure (startTime + exposureTime + readoutTime)
 *
 * \pre Sub-pixel, pixel, bias register, and smearing map filled with values from previous exposure.
 *
 * \post Pixel unit in the pixel, bias register, and smearing maps: [ADU]
 */

double DetectorWithMappedPSF::takeExposure(int exposureNr, double startTime, double exposureTime)
{
    // Advance the internal clock until the given start time

    internalTime = startTime;

    // Integration of point sources and background, taking into account jitter + drift.

    Log.info("Detector: Integrating light for exposure " + to_string(exposureNr) + " with exposure time = " + to_string(exposureTime));

    integrateLight(exposureNr, startTime, exposureTime);

    // Include noise effects like readout noise, photon noise, full well saturation, etc.
    // Note: readOut() needs the exposure time to compute the open shutter smearing.

    Log.info("Detector: Adding noise effects to exposure " + to_string(exposureNr));

    readOut(exposureTime);

    // Write the CCD subfield, the bias map, and the smearing map to the HDF5 file

    Log.debug("Detector: Writing PixelMap, smearing map, and bias map #" + to_string(exposureNr) + " to HDF5 file.");

    writePixelMapsToHDF5(exposureNr);

    // If required, also write the subpixel image to the HDF5 file

    if (writeSubPixelImagesToHDF5)
    {
        Log.debug("Detector: Writing SubPixelMap " + to_string(exposureNr) + " to HDF5 file.");
        writeSubPixelMapToHDF5(exposureNr);
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
 * \param exposureNr:   Sequential number of the exposure
 * \param startTime:    Starting time of the exposure for which jitter must be applied [s].
 * \param exposureTime: Duration of the exposure [s].
 *
 * \pre Sub-pixel, pixel, bias register, and smearing map filled with values from previous exposure.
 *
 * \post Pixel unit of the sub-pixel map: [photons].
 * \post Pixel, bias register, and smearing map filled with zeroes.
 */

void DetectorWithMappedPSF::integrateLight(int exposureNr, double startTime, double exposureTime)
{

    // Reset the sub-field (i.e. get rid of the previous exposure, by zeroing the entire sub-field)

    Log.debug("Detector: resetting subfield array for new exposure.");

    reset();

    // Integration (incl. jitter): point sources + background

    camera.exposeDetector(*this, startTime, exposureTime);

    // Convolve with the point spread function

    convolveWithPsf();

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
    // This takes into account the QE, vignetting, polarisation, and particulate & molecular contamination.
    // PixelMap units change from [photons] to [electrons] 

    applyThroughputEfficiency();

    // Add dark current

    if(includeDarkSignal)
    {
    		Log.debug("Detector: adding dark current");

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
 * \param xFP   X-coordinate of the sub-pixel in the focal plane in the FP reference frame [mm].
 * \param yFP   Y-coordinate of the sub-pixel in the focal plane in the FP reference frame [mm].
 * \param flux  Flux to add to the sub-pixel map [photons].
 *
 * \return           (isInSubfield, row, col) 
 *                   isInSubfield: True if (xFP, yFP) are on the subfield, false otherwise.
 *                   row:          subfield (not CCD) row number of the pixel to which the flux was added
 *                   col:          subfield (not CCD) column number of the pixel to which the flux was added  
 */

tuple<bool, double, double> DetectorWithMappedPSF::addFlux(double xFP, double yFP, double flux)
{
    // Convert from FP coordinates to CCD pixel coordinates

    double pixRow, pixColumn;
    tie(pixRow, pixColumn) = focalPlaneToPixelCoordinates(xFP, yFP);

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

bool DetectorWithMappedPSF::isInSubPixelMap(double row, double column)
{
    return (column >= 0) && (row >= 0) && (column < numColumnsSubPixelMap) && (row < numRowsSubPixelMap);
}











/**
 * \brief: Add the given flux value to (all sub-pixels of) the sub-pixel map.
 *
 * \param flux: Flux to add to the sub-pixel map [photons/pixel].
 *
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
 * \pre Flatfield map at sub-pixel level, excl. edge pixels.
 * \pre Pixel, bias register, and smearing maps filled with zeroes.
 *
 * \post Pixel value in the sub-pixel map: [photons].
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
 * \brief: Rebin the sub-pixel map to pixel level and crop the edge pixels.
 *
 * \pre Unit of the pixel value in the sub-pixel map: [photons].
 * \pre Pixel, bias register, and smearing map filled with zeroes.
 *
 * \post Unit of pixel values in the sub-pixel map: [photons].
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

void DetectorWithMappedPSF::setPsfForSubfield()
{
    // There is one PSF for the entire subfield, which we take the one of the center 
    // of the subfield.

    double xFPmm, yFPmm;
    tie(xFPmm, yFPmm) = getFocalPlaneCoordinatesOfSubfieldCenter();

    // Get the 'user specified' angular distance to the optical axis from the psf.
    // If the user didn't specify an angular distance, calculate it from the given
    // focal plane coordinates.

    double radius = psf->getRequestedDistanceToOpticalAxis();

    if (radius < 0.0)
    {
        radius = camera.getGnomonicRadialDistanceFromOpticalAxis(xFPmm, yFPmm);
    }

    psf->select(radius);


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
    psf->rotate(angle);

    // Rebin the psfMap to the number of sub-pixels per pixel used for the Detector

    psfMap = psf->rebinToSubPixels(numSubPixelsPerPixel);

    // Allow the convolver to precompute some stuff given the PSF, so that it doesn't
    // need to be recomputed every convolution.

    convolver.initialise(numRowsSubPixelMap, numColumnsSubPixelMap, psfMap);

}











/**
 * \brief: Convolve the sub-pixel map with the PSF, keeping the same dimensions.
 *
 * \param psf: PSF.
 */

void DetectorWithMappedPSF::convolveWithPsf()
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




