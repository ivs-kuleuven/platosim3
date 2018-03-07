#include "DetectorWithAnalyticNonGaussianPSF.h"





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
 * The flatfieldMap, throughputMap and cteMap are filled at pixel level.
 * 
 * \param configParam    Configuration parameters for the detector.
 * \param hdf5file       HFD5 file to write the detector images to.
 * \param camera         Camera to which to attach the detector.
 */

DetectorWithAnalyticNonGaussianPSF::DetectorWithAnalyticNonGaussianPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator)
: Detector(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator), sigma(nullptr)
{
    // Parse the parameters from the configuration file.

    configure(configParam);

    // Allocate memory for the flatfield map

    flatfieldMap.ones(numRowsPixelMap, numColumnsPixelMap);

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
}








/**
 * Destructor.
 *
 */
DetectorWithAnalyticNonGaussianPSF::~DetectorWithAnalyticNonGaussianPSF()
{
    flushOutput();
    delete sigma;
}










/**
 * \brief Configure the DetectorWithAnalyticNonGaussianPSF object using the ConfigurationParameters
 * 
 * \param configParam: the configuration parameters 
 **/

void DetectorWithAnalyticNonGaussianPSF::configure(ConfigurationParameters &configParam)
{
    flatfieldNoiseAmplitude   = configParam.getDouble("CCD/FlatfieldPtPNoise");
    includeFlatfield          = configParam.getBoolean("CCD/IncludeFlatfield");
    flatfieldSeed             = configParam.getLong("RandomSeeds/FlatFieldSeed");


    string filename = configParam.getAbsoluteFilename("PSF/AnalyticNonGaussian/ParameterFileName");

    ifstream file(filename);
    if (!file)
        return;

    params.clear();
    string line;
    while (getline(file, line)) 
    {
        if (line == "" || line.find("#") == 0)
            continue;

        istringstream strs(line);
        vector<double> vals((istream_iterator<double>(strs)), istream_iterator<double>());

        if ((params.size() > 0 && params[0].size() == vals.size()) || (params.size() == 0 && vals.size() > 0))
        {
            params.emplace_back(vals);
        }
    }
    
    // The sigma of the PSF can either be a fixed value, or given by a time series in a file
    
    string sigmaPSFSource = configParam.getString("PSF/AnalyticNonGaussian/Sigma/Source");
    if (sigmaPSFSource == "ConstantValue")
    {
        double sigmaPSFValue = configParam.getDouble("PSF/AnalyticNonGaussian/Sigma/ConstantValue");     // [pix]
        sigma = new Parameter<double>(sigmaPSFValue);
    
        Log.info("DetectorWithAnalyticNonGaussianPSF: Using a constant sigma: " + to_string(sigmaPSFValue) + " pix");
    }
    else if (sigmaPSFSource == "FromFile")
    {
        string sigmaPSFInputFile = configParam.getAbsoluteFilename("PSF/AnalyticNonGaussian/Sigma/FromFile");
        sigma = new Parameter<double>(sigmaPSFInputFile, 1);                                            // [pix]
    
        Log.info("DetectorWithAnalyticNonGaussianPSF: Reading sigma PSF from " + sigmaPSFInputFile);
    }
}












/**
 * \brief Update the time dependent parameters of the Detector to their 
 *        value at the given time point
 *
 * \param time: current time
 *
 * \return 
 */

void DetectorWithAnalyticNonGaussianPSF::updateParameters(double time)
{
    sigma->updateValue(time);
}











/**
 * \brief Interpolate and rotate PSF parameters and sum up all parts to calculate the intergal of the analytic PSF.
 * 
 * \param psf:   container to hold the result of the integration
 * \param x:     x position of the PSF
 * \param y:     y position of the PSF
 * \param r:     radial distance of the PSF to the optical axis
 * \param p:     azimuth angle of the PSF
 * \param scale: scale factor to resize the PSF
 **/

void DetectorWithAnalyticNonGaussianPSF::integrateAnalyticPSF(IntegralOfAnalyticSignalResponse& psf, double x, double y, double r, double p, double scale) 
{
    double ox = x - floor(x);
    double oy = y - floor(y);
    double s = (*sigma)() * scale;
    if (params.size() > 0 && params[0].size() > 6) 
    {
        r /= 1.4;
        unsigned c1 = min(params[0].size() / 7 - 1, (size_t)r) * 7;
        unsigned c2 = min(params[0].size() / 7 - 1, (size_t)r + 1) * 7;
        double w = r - (unsigned)r;
        w = 3. * w * w - 2. * w * w * w;

        for (auto i = params.cbegin(); i != params.cend(); i++) 
        {
            double pr = s * ((1. - w) * (*i)[c1] + w * (*i)[c2]);
            double pp = (1. - w) * (*i)[c1 + 1] + w * (*i)[c2 + 1];
            double h = (1. - w) * (*i)[c1 + 2] + w * (*i)[c2 + 2];
            double b = s * ((1. - w) * (*i)[c1 + 3] + w * (*i)[c2 + 3]);
            double r = s * ((1. - w) * (*i)[c1 + 4] + w * (*i)[c2 + 4]);
            double m = (1. - w) * (*i)[c1 + 5] + w * (*i)[c2 + 5];
            double a = (1. - w) * (*i)[c1 + 6] + w * (*i)[c2 + 6];

            psf.addPart(ox + pr * cos(p + pp), oy + pr * sin(p + pp), h, b, r, m, p + a);
            psf.addPart(ox + pr * cos(p - pp), oy + pr * sin(p - pp), h, b, r, m, p - a);
        }
    } 
    else 
    {
        psf.addPart(ox, oy, 1., s);
    }
}











 /**
 * \brief: Generate the (random) flatfield variations.  This map is generated
 *         at pixel level but without the edge pixels.
 *
 * https://github.com/python-acoustics/python-acoustics/blob/master/acoustics/generator.py#L108
 */

void DetectorWithAnalyticNonGaussianPSF::generateFlatfieldMap()
{

    Log.info("Detector: generating flatfield map.");

    // Random number generation

    mt19937 flatfieldGenerator(flatfieldSeed);
    normal_distribution<double> flatfieldDistribution(0.0, 1.0);

    // Double the dimensions (this is necessary because of the behaviour of the Fourier transforms)

    int Nrows = 2 * numRowsPixelMap;
    int Ncolumns = 2 * numColumnsPixelMap;

    arma::cx_fmat evenMap = arma::cx_fmat(Nrows, Ncolumns);

    for(unsigned int row = 0; row < Nrows; row++)
    {
        for(unsigned int column = 0; column < Ncolumns; column++)
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

    flatfieldMap(arma::span::all, arma::span::all) = realMap(arma::span(0, Nrows / 2 - 1), arma::span(0, Ncolumns / 2 - 1));

    // Normalise

    float minPinkNoise = flatfieldMap.min();
    float maxPinkNoise = flatfieldMap.max();

    flatfieldMap -= minPinkNoise;
    flatfieldMap /= (maxPinkNoise - minPinkNoise); // [0, 1]
    flatfieldMap *= flatfieldNoiseAmplitude;    // [0, flatfialdNoiseAmplitude]
    flatfieldMap += (1.0 - flatfieldNoiseAmplitude);

    // Write the result to the HDF5 output file

    Log.debug("Detector: writing PRNU to HDF5");

    hdf5File.writeArray("/Flatfield", "PRNU", flatfieldMap);
}













/**
 * \brief: Zeroes the pixel, bias register, and the smearing maps.
 *
 * \pre pixel, bias register, and smearing maps filled with values from previous exposure.
 *
 * \post pixel, bias register, and smearing maps filled with zeroes.
 */

void DetectorWithAnalyticNonGaussianPSF::reset()
{
    pixelMap.zeros();
    biasMap.zeros();
    smearingMap.zeros();
}










/**
 * \brief: Take an exposure with the detector starting at the given time.
 *         The light is integrated during the given exposure time, during which 
 *         the detector experiences the effects of jitter and thermo-elastic telescope 
 *         drift. The background is assumed uniform for the whole subfield.
 *         Afterwards, the collected light is read out, and various noise effects are added.
 *
 * \param exposureNr:   Sequential number of the exposure
 * \param startTime:    Starting time of the exposure [s].
 * \param exposureTime: Duration of the exposure [s].
 * 
 * \return endTime:     Time after the exposure (startTime + exposureTime + readoutTime)
 *
 * \pre Pixel, bias register, and smearing map filled with values from previous exposure.
 *
 * \post Pixel unit in the pixel, bias register, and smearing maps: [ADU]
 */

double DetectorWithAnalyticNonGaussianPSF::takeExposure(int exposureNr, double startTime, double exposureTime)
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
 * \param startTime:    Starting time of the exposure [s].
 * \param exposureTime: Duration of the exposure [s].
 *
 * \pre Pixel, bias register, and smearing map filled with values from previous exposure.
 *
 * \post Pixel unit of the sub-pixel map: [photons].
 * \post Pixel, bias register, and smearing map filled with zeroes.
 */

void DetectorWithAnalyticNonGaussianPSF::integrateLight(int exposureNr, double startTime, double exposureTime)
{

    // Reset the subfield (i.e. get rid of the previous exposure, by zeroing the entire sub-field)

    Log.debug("Detector: resetting subfield array for new exposure.");

    reset();

    // Integration (incl. jitter): point sources + background

    camera.exposeDetector(*this, startTime, exposureTime);

    // Apply flatfield (at pixel level)

    if (includeFlatfield)
    {
        Log.debug("Detector: applying Flatfield.");

        applyFlatfield();
    }
    else
    {
        Log.debug("Detector: no flatfield applied.");
    }

    // Apply throughput efficiency on the pixel map
    // This takes into account the QE, vignetting, polarisation, and particulate & molecular contamination.
    // PixelMap units change from [photons] to [electrons] 

    applyThroughputEfficiency();

    // Brighter-Fatter effect

    if(includeBFE)
    {
    		Log.debug("Detector: adding Brighter-Fatter effect");

       	applyBFE();
    }

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
 * \brief: Add the PSF of the star with given focal plane coordinates and flux level to the pixel map.
 *         Return the pixel coordinates of the barycenter of the PSF. As PSF we use an analytic non-Gaussian 
 *         function.
 *         
 * \param xFP   X-coordinate of the (fractional) pixel in the focal plane in the FP reference frame [mm].
 * \param yFP   Y-coordinate of the (fractional) pixel in the focal plane in the FP reference frame [mm].
 * \param flux  Flux to add to the pixel map [photons].
 *
 * \return           (isInSubfield, row, col) 
 *                   isInSubfield: True if (xFP, yFP) are on the subfield, false otherwise.
 *                   row:          subfield (not CCD) row number of the barycenter of the PSF.
 *                   col:          subfield (not CCD) column number of the barycenter of the PSF.
 *                   
 * \note In the code below, we use pixelMap(i,j) that involves a boundary check. To increase speed,
 *       switch to pixelMap.at(i,j) that does not involve a boundary check. 
 */

tuple<bool, double, double> DetectorWithAnalyticNonGaussianPSF::addFlux(double xFP, double yFP, double flux)
{
    // Convert the FP coordinates of the PSF barycenter to (real-valued) CCD pixel coordinates
    // Then convert from CCD pixel coordinates to subfield pixel coordinates.

    double row0, column0;
    tie(row0, column0) = focalPlaneToPixelCoordinates(xFP, yFP);
    row0 -= subFieldZeroPointRow;
    column0 -= subFieldZeroPointColumn;

    int size = 2 * ((int)(8. * (*sigma)()) + 1) + 1;;
    int sx = (int)floor(column0 - (size - 1.) / 2.);
    int sy = (int)floor(row0 - (size - 1.) / 2.);

    if (sx + size <= 0 || sx >= (int)numColumnsPixelMap || sy + size <= 0 || sy >= (int)numRowsPixelMap)
        return make_tuple(false, row0, column0);

    IntegralOfAnalyticSignalResponse psf(size);
    double r = rad2deg(camera.getGnomonicRadialDistanceFromOpticalAxis(xFP, yFP));
    double p = atan2(yFP, xFP);

    double ccdOrientation = getOrientationAngle();
    p -= ccdOrientation;

    integrateAnalyticPSF(psf, column0, row0, r, p);

    for (int y = max(0, sy); y < min((int)numRowsPixelMap, sy + size); y++)
        for (int x = max(0, sx); x < min((int)numColumnsPixelMap, sx + size); x++)
            pixelMap.at(y, x) += psf(x - sx, y - sy) * flux;

    return make_tuple(true, row0, column0);
}












/**
 * \brief   Check whether the given (row, column) indices are within the array range of the pixel map.
 *
 * \details  The input parameters row & column come from a coordinate transformation
 *           in the focal plane, and as a result are not necessarily integers. For this 
 *           function it's not necessary to round them to the nearest integer. 
 *
 * \param  row:    Row index. NOT a coordinate in the CCD frame, but in the subfield frame.    [pixel].
 * \param  column: Column index. NOT a coordinate in the CCD frame, but in the subfield frame. [pixel].
 *
 * \return  True if the given (row, column) coordinates are in the pixel map; false otherwise.
 */

bool DetectorWithAnalyticNonGaussianPSF::isInPixelMap(double row, double column)
{
    return (column >= 0) && (row >= 0) && (column < numColumnsPixelMap) && (row < numRowsPixelMap);
}











/**
 * \brief: Add the given flux value to (all sub-pixels of) the sub-pixel map.
 *
 * \param flux: Flux to add to the sub-pixel map [photons/pixel].
 *
 */

void DetectorWithAnalyticNonGaussianPSF::addFlux(double flux)
{
    pixelMap += flux;
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

void DetectorWithAnalyticNonGaussianPSF::applyFlatfield()
{
    const unsigned int beginRow = numEdgePixels;
    const unsigned int beginCol = numEdgePixels;
    const unsigned int endRow = numRowsPixelMap - numEdgePixels - 1;
    const unsigned int endCol = numColumnsPixelMap - numEdgePixels - 1;

    pixelMap.submat(beginRow, beginCol, endRow, endCol) = pixelMap.submat(beginRow, beginCol, endRow, endCol) % flatfieldMap;
}

