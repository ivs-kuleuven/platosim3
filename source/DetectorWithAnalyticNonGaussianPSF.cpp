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
 * \param readoutTimeBeforeNextExposure Duration of the readout that takes place before the next exposure can start.
 */

DetectorWithAnalyticNonGaussianPSF::DetectorWithAnalyticNonGaussianPSF(ConfigurationParameters &configParam,
								       HDF5File &hdf5file,
								       Camera &camera,
								       TemperatureGenerator &feeTemperatureGenerator,
								       TemperatureGenerator &detectorTemperatureGenerator,
								       double readoutTimeBeforeNextExposure,
								       double readoutTimeDuringNextExposure)
  : Detector(configParam,
	     hdf5file,
	     camera,
	     feeTemperatureGenerator,
	     detectorTemperatureGenerator,
	     readoutTimeBeforeNextExposure,
	     readoutTimeDuringNextExposure),
    sigma(nullptr)
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
    numExposures        = configParam.getUnsignedInteger("ObservingParameters/NumExposures");
    beginExposureNr     = configParam.getUnsignedInteger("ObservingParameters/BeginExposureNr");
    cycleTime           = configParam.getDouble("ObservingParameters/CycleTime");
    flatfieldNoiseRMS   = configParam.getDouble("CCD/FlatfieldNoiseRMS");
    includeFlatfield    = configParam.getBoolean("CCD/IncludeFlatfield");
    flatfieldSeed       = configParam.getLong("RandomSeeds/FlatFieldSeed");

    // Read and configure the parameters used to calculate the PSF

    string filename = configParam.getAbsoluteFilename("PSF/AnalyticNonGaussian/ParameterFileName");

    ifstream file(filename);
    if (!file)
    {
        Log.error("DetectorWithAnalyticNonGaussianPSF::configure(): Parameter file doesn't exist or is not readable: "  + filename);
        throw ConfigurationException("DetectorWithAnalyticNonGaussianPSF: wrong parameter filename in configuration file");
    }

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

    // The parameters for the charge diffusion

    includeChargeDiffusion = configParam.getBoolean("PSF/AnalyticNonGaussian/IncludeChargeDiffusion");
    chargeDiffusionStrength = configParam.getDouble("PSF/AnalyticNonGaussian/ChargeDiffusionStrength");

    Log.info("DetectorWithAnalyticNonGaussianPSF: sigma of charge diffusion: " + to_string(chargeDiffusionStrength) + " pix");

    // The sigma of the PSF can either be a fixed value, or given by a time series in a file

    string sigmaPSFSource = configParam.getString("PSF/AnalyticNonGaussian/Sigma/Source");
    if (sigmaPSFSource == "ConstantValue")
    {
        double sigmaPSFValue = configParam.getDouble("PSF/AnalyticNonGaussian/Sigma/ConstantValue");     // [pix]
        sigma = new Parameter<double>(sigmaPSFValue);

        Log.info("DetectorWithAnalyticNonGaussianPSF: Using a constant PSF sigma: " + to_string(sigmaPSFValue) + " pix");
    }
    else if (sigmaPSFSource == "FromFile")
    {
        string sigmaPSFInputFile = configParam.getAbsoluteFilename("PSF/AnalyticNonGaussian/Sigma/FromFile");
        sigma = new Parameter<double>(sigmaPSFInputFile, 1);                                            // [pix]

        Log.info("DetectorWithAnalyticNonGaussianPSF: Reading sigma PSF from " + sigmaPSFInputFile);
    }

    // The configuration for the on-the-fly photometry

    includePhotometry = configParam.getBoolean("Photometry/IncludePhotometry");

    if (includePhotometry)
    {
        photometry = Photometry(configParam, hdf5file, camera);
    }

    // The configuration for the HDF5 contents

    writeFlatfieldMap      = configParam.getBoolean("ControlHDF5Content/WriteFlatfieldMap");
    writeHighResolutionPSF = configParam.getBoolean("ControlHDF5Content/WriteHighResolutionPSF");
    writeDiffusedPSF       = configParam.getBoolean("ControlHDF5Content/WriteDiffusedPSF");

} // end configure()





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
    if (ccdPosition != "Custom")
        ccdPositions->updateValue(time);

    sigma->updateValue(time);
}





/**
 * \brief Interpolate and rotate PSF parameters and sum up all parts to calculate the integral of the analytic PSF.
 *
 * \param psf:        container to hold the result of the integration
 * \param x:          x position of the PSF
 * \param y:          y position of the PSF
 * \param r:          radial distance of the PSF to the optical axis
 * \param p:          azimuth angle of the PSF
 * \param Nsubpixels: number of subpixels per pixel (e.g. 128, set to 1 for no subpixels)
 **/

void DetectorWithAnalyticNonGaussianPSF::integrateAnalyticPSF(IntegralOfAnalyticSignalResponse& psf, double x, double y, double q, double p, double rotpsf, int Nsubpixels)
{
    double ox = x - floor(x);
    double oy = y - floor(y);
    double s = (*sigma)() * Nsubpixels;
    if (params.size() == 6 && params[0].size() > 6)
    {
        q /= 1.4;
        p -= rotpsf;
        unsigned c1 = min(params[0].size() / 7 - 1, (size_t)q) * 7;
        unsigned c2 = min(params[0].size() / 7 - 1, (size_t)q + 1) * 7;
        double w = q - (unsigned)q;
        w = 3. * w * w - 2. * w * w * w;

        for (auto i = params.cbegin(); i != params.cend(); i++)
        {
            double pr = s * ((1. - w) * (*i)[c1] + w * (*i)[c2]);
            double pp = (1. - w) * (*i)[c1 + 1] + w * (*i)[c2 + 1];
            double h  = (1. - w) * (*i)[c1 + 2] + w * (*i)[c2 + 2];
            double b  = s * ((1. - w) * (*i)[c1 + 3] + w * (*i)[c2 + 3]);
            double rr = s * ((1. - w) * (*i)[c1 + 4] + w * (*i)[c2 + 4]);
            double m  = (1. - w) * (*i)[c1 + 5] + w * (*i)[c2 + 5];
            double a  = (1. - w) * (*i)[c1 + 6] + w * (*i)[c2 + 6];

            psf.addPart(ox + pr * cos(p + pp), oy + pr * sin(p + pp), h, b, rr, m, p + a);
            psf.addPart(ox + pr * cos(p - pp), oy + pr * sin(p - pp), h, b, rr, m, p - a);
        }
    }
    else if (params.size() > 0 && params[0].size() > 6)
    {
        unsigned S = params[0].size() / 7;
        unsigned N = params.size() / S;
        double c = atan(tan(q / 180. * M_PI) * cos(p)) / M_PI * 180. * S / 36. + (S - 1.) / 2.;
        double r = atan(tan(q / 180. * M_PI) * sin(p)) / M_PI * 180. * S / 36. + (S - 1.) / 2.;
        double t = p - rotpsf;
        double co = s * cos(t), si = s * sin(t);
        double wc = c - floor(c);
        double wr = r - floor(r);
        unsigned c1 = (unsigned)min(max(c, 0.), S - 1.) * 7;
        unsigned c2 = (unsigned)min(max(c + 1., 0.), S - 1.) * 7;
        unsigned r1 = (unsigned)min(max(r, 0.), S - 1.) * N;
        unsigned r2 = (unsigned)min(max(r + 1., 0.), S - 1.) * N;
        double d[7];

        for (unsigned n = 0; n < N; n++)
        {
            for (unsigned i = 0; i < 7; i++)
            {
                double d1 = (1. - wc) * params[r1 + n][c1 + i] + wc * params[r1 + n][c2 + i];
                double d2 = (1. - wc) * params[r2 + n][c1 + i] + wc * params[r2 + n][c2 + i];
                d[i] = (1. - wr) * d1 + wr * d2;
            }
            psf.addPart(ox + d[0] * co - d[1] * si, oy + d[0] * si + d[1] * co, d[2], d[3] * s, d[4] * s, d[5], d[6] + t);
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
        Log.debug("Detector: writing PRNU to HDF5");
        hdf5File.writeArray("/Flatfield", "PRNU", flatfieldMap);
    }
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

    // Clear all arrays

    Log.debug("Detector: resetting subfield array for new exposure.");
    reset();

    // Integration of point sources and background, taking into account jitter + drift.
    
    if (!constantSkyBackground && (exposureNr == beginExposureNr))
    {
        fillBackgroundMap(camera, startTime, exposureTime);
    }

    Log.info("DetectorWithAnalyticNonGaussianPSF: Integrating light for exposure " + to_string(exposureNr) + " with exposure time = " + to_string(exposureTime));

    integrateLight(exposureNr, startTime, exposureTime);

    // Include noise effects like readout noise, photon noise, full well saturation, etc.
    // Note: readOut() needs the exposure time to compute the open shutter smearing.

    Log.info("DetectorWithAnalyticNonGaussianPSF: Adding noise effects to exposure " + to_string(exposureNr));

    if (exposureNr == beginExposureNr) {
      if (includeCTIeffects &&
          (CTImodel == "Short2013" || CTImodel == "Short2013FromFile"))
      {
          setInitialNumberOfOccupiedTraps(numberOfOccupiedTrapsPixelMap);
      }
    }

    readOut(exposureTime);

    // If photometric extraction was asked, apply it now

    if (includePhotometry)
    {
        Log.info("DetectorWithAnalyticNonGaussianPSF: applying photometric extraction to exposure " + to_string(exposureNr));
        photometry.applyPhotometry(exposureNr);
    }

    // Write the CCD subfield, the bias map, and the smearing map to the HDF5 file

    Log.debug("DetectorWithAnalyticNonGaussianPSF: Writing PixelMap, smearing map, and bias map #" + to_string(exposureNr) + " to HDF5 file.");

    writePixelMapsToHDF5(exposureNr);

    // Write the cosmic hits to the HDF5 file

    Log.debug("DetectorWithAnalyticNonGaussianPSF: Writing cosmics of the PixelMap, smearing map, bias map #" + to_string(exposureNr) + " to HDF5 file.");

    if (writeCosmics)
    {
            if (groupByExposure){writeCosmicHitsToHDF5WhenGroupByExposure(exposureNr);}
            else{writeCosmicHitsToHDF5WithoutGroupByExposure(exposureNr);}
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

    // Integration (incl. jitter): point sources + background

    camera.exposeDetectorWithStars(*this, startTime, exposureTime, readoutTimeBeforeNextExposure);
    if (constantSkyBackground)
    {
        camera.exposeDetectorWithSkyBackground(*this, startTime, exposureTime, readoutTimeBeforeNextExposure);
    }
    else
    {
        addBackgroundMapToPixelMap(camera, startTime);
    }
    // Apply throughput efficiency on the pixel map.
    // This takes into account the QE, vignetting, polarisation, and particulate & molecular contamination.
    // PixelMap units change from [photons] to [electrons]

    applyThroughputEfficiency();

    // Apply the charge injection which will mitigate the CTI. The injection happens in electrons,
    // so the throughput efficiency should already have been applied. The injected charges do feel the PRNU,
    // so applying the flatfied should happen afterwards.

    if (includeChargeInjection)
    {
        Log.debug("Detector: applying charge injection");
        applyChargeInjection();
    }

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
 * \brief: Add the PSF of the star with given focal plane coordinates and flux level to the given map.
 *         As PSF we use an analytic non-Gaussian function. This function gets called in the addFlux()
 *         method to add the flux to the pixelMap, and in the applyPhotometry() function. This method
 *         does not exist in any other child class of the detector class.
 *
 * \param map      matrix with the same dimensions as pixelMap
 * \param row0     real-valued subfield row index of the star position               [pix]
 * \param column0  real-values subfield column index of the star position            [pix]
 * \param r        angular distance of the star position to the optical axis         [rad]
 * \param p        azimuthal angle in the focal plane of the star position           [rad]
 * \param flux     flux to add to the  map                                           [photons]
 *
 * \return         isInSubfield: true if the entire PSF is on the map, false otherwise.
 */

bool DetectorWithAnalyticNonGaussianPSF::addFluxToMap(arma::Mat<float>& map, double row0, double column0, double r, double p, double flux)
{
    double s = (*sigma)();
    double diffusionKernelWidth = 0.;

    if (includeChargeDiffusion)
    {
        diffusionKernelWidth = chargeDiffusionStrength;
        s = sqrt(s * s + diffusionKernelWidth * diffusionKernelWidth);
    }

    int size = 2 * (int)(8. * s + 1) + 1;
    int sx = (int)floor(column0 - (size - 1.) / 2.);
    int sy = (int)floor(row0 - (size - 1.) / 2.);

    if (sx + size <= 0 || sx >= (int)map.n_cols || sy + size <= 0 || sy >= (int)map.n_rows)
        return false;

    // Construct the PSF around the central pixel coordinates

    IntegralOfAnalyticSignalResponse psf(size, diffusionKernelWidth);

    integrateAnalyticPSF(psf, column0, row0, r, p, rotationAnglePsf);

    for (int y = max(0, sy); y < min((int)map.n_rows, sy + size); y++)
        for (int x = max(0, sx); x < min((int)map.n_cols, sx + size); x++)
            map.at(y, x) += psf(x - sx, y - sy) * flux;

    return true;
}





/**
 * \brief: Add the PSF of the star with given focal plane coordinates and flux level to the pixel map.
 *         Return the pixel coordinates of the barycenter of the PSF. As PSF we use an analytic non-Gaussian
 *         function. The flux gets added to the pixelMap using the fuction addFluxToMap().
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
 */

tuple<bool, double, double> DetectorWithAnalyticNonGaussianPSF::addFlux(double xFP, double yFP, double flux)
{
    // Convert the FP coordinates of the PSF barycenter to (real-valued) CCD pixel coordinates
    // Then convert from CCD pixel coordinates to subfield pixel coordinates.

    double row0, column0;
    tie(row0, column0) = focalPlaneToPixelCoordinates(xFP, yFP);
    row0 -= subFieldZeroPointRow;
    column0 -= subFieldZeroPointColumn;

    // Check if the star falls in the subfield. If not, don't add any flux, but simply return.

    if (!isInPixelMap(row0, column0))
    {
        return make_tuple(false, row0, column0);
    }

    // Get the polar coordinates of the star in the focal plane. This determines the shape of the PSF.

    double r = rad2deg(camera.getGnomonicRadialDistanceFromOpticalAxis(xFP, yFP));
    double p = atan2(yFP, xFP);

    // Add the PSF to the pixel map. If the PSF does not entirely fall on the pixel map, then success == false.

    bool success = addFluxToMap(pixelMap, row0, column0, r, p, flux);

    return  make_tuple(success, row0, column0);
}





/**
 * \brief Insert the extended ghost with the given radius and flux at the given focal-plane position.
 *
 * Note that the extended source will not be convolved with the PSF, for practical reasons (but since the
 * extended ghosts are so large, the influence of the PSF is negligible).
 *
 * \param x0: Focal-plane x-coordinate of the centre of the extended ghost [mm].
 * \param y0: Focal-plane y-coordinate of the centre of the extended ghost [mm].
 * \param radius: Radius of the extended ghost [mm].
 * \param flux: Flux of the extended ghost [photons].
 *
 * \return: Whether or not the extended source falls (at least partially) on the sub-field, and the
 *          (row, column) coordinates of the centre of the extended ghost in the pixel map.
 */
tuple<bool, double, double> DetectorWithAnalyticNonGaussianPSF::addExtendedGhost(double x0, double y0, double radius, double flux)
{
    // Calculate the number of pixels in the extended ghost

    double radiusPixels = radius * 1000 / pixelSize;    // Radius [pixels]
    double radiusPixelsSquared = pow(radiusPixels, 2);  // Squared radius [pixels^2]

    double numPixels = PI * pow(radiusPixels, 2);       // Area of the extended ghost [pixels]
    double fluxPerPixel = flux / numPixels;             // Flux [photons / pixel]

    // Calculate the (row, column) coordinates of the centre of the extended source in the pixel map

    double row0, column0;
    tie(row0, column0) = focalPlaneToPixelCoordinates(x0, y0);
    row0 -= subFieldZeroPointRow;
    column0 -= subFieldZeroPointColumn;

    bool ghostInPixelMap = false;

    // Try to add flux to all pixels covered by the extended ghosts

    for(int row = row0 - radiusPixels; row <= row0 + radiusPixels; row++)
    {
        for(int column = column0 - radiusPixels; column <= column0 + radiusPixels; column++)
        {
            if (isInPixelMap(row, column) && pow(column - column0, 2) + pow(row - row0, 2) <= radiusPixelsSquared)
            {
                ghostInPixelMap = true;
                pixelMap(row, column) += fluxPerPixel;
            }
        }
    }

    return  make_tuple(ghostInPixelMap, row0, column0);
}





/**
 * \brief: Create a high-resolution map of the PSF in the center of the subfield
 *
 * \param map         The PSF will be written in this 2D array, the size need not be allocated.
 * \param Npixels     Size of the high-res map in pixels
 * \param Nsubpixels  (sqrt of the) number of subpixels per pixel
 *
 * \return           'map' will be modified to contain the high-resolution PSF
 *                   Its size will be  (Npixels * Nsubpixels) x  (Npixels * Nsubpixels)
 */

void DetectorWithAnalyticNonGaussianPSF::makeHighResolutionPSF(arma::Mat<float> &highResMap, bool includeDiffusion, int Npixels, int Nsubpixels)
{
    // Put the PSF right in the middle of the (high-res) (sub)pixel map

    double row0    = Npixels / 2.0;
    double column0 = Npixels / 2.0;

    // The high-res (sub)pixel map will be placed in the middle of the regular subfield.
    // Derive the focal plane coordinates of the middle of the subfield. These are needed
    // to later on derive the angular distance from the optical axis.

    double middleRowSubfield = subFieldZeroPointRow    + numRowsPixelMap    / 2.0;
    double middleColSubfield = subFieldZeroPointColumn + numColumnsPixelMap / 2.0;
    double xFP, yFP;
    tie(xFP, yFP) = pixelToFocalPlaneCoordinates(middleRowSubfield, middleColSubfield);

    // Convolve with Gaussian diffusion kernel is requested

    double diffusionKernelWidth = 0.;

    if (includeDiffusion)
    {
        diffusionKernelWidth = chargeDiffusionStrength;
    }

    int size = Npixels * Nsubpixels;
    highResMap.set_size(size, size);
    highResMap.fill(0.0);

    int sx = (int)floor(column0 * Nsubpixels - (size - 1.) / 2.);
    int sy = (int)floor(row0    * Nsubpixels - (size - 1.) / 2.);

    // Construct the PSF around the central pixel coordinates

    IntegralOfAnalyticSignalResponse psf(size, diffusionKernelWidth*Nsubpixels);
    double r = rad2deg(camera.getGnomonicRadialDistanceFromOpticalAxis(xFP, yFP));
    double p = atan2(yFP, xFP);

    integrateAnalyticPSF(psf, column0, row0, r, p, rotationAnglePsf, Nsubpixels);

    for (int y = max(0, sy); y < min((int)Npixels*Nsubpixels, sy + size); y++)
        for (int x = max(0, sx); x < min((int)Npixels*Nsubpixels, sx + size); x++)
            highResMap.at(y, x) += psf(x - sx, y - sy) * 1000.0;

    // Normalize the PSF

    highResMap /= arma::accu(highResMap);
}





/**
 * \brief: Add the given flux value to (all sub-pixels of) the sub-pixel map.
 *
 * \param flux: Flux to add to the sub-pixel map [photons/pixel].
 *
 */

void DetectorWithAnalyticNonGaussianPSF::addFlux(double flux)
{
    bool subFieldIsBlockedOff = (coveredBottom + coveredTop >= numRowsPixelMap) || (coveredLeft + coveredRight >= numColumnsPixelMap);
    if (!subFieldIsBlockedOff)
    {
        pixelMap.submat(coveredBottom, coveredLeft,
                      numRowsPixelMap - coveredTop - 1,
                      numColumnsPixelMap - coveredRight - 1) += flux;
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
    const unsigned int endRow   = numRowsPixelMap    - numEdgePixels - 1;
    const unsigned int endCol   = numColumnsPixelMap - numEdgePixels - 1;

    pixelMap.submat(beginRow, beginCol, endRow, endCol) = pixelMap.submat(beginRow, beginCol, endRow, endCol) % flatfieldMap;
}





/**
 *  \brief Before destroying this object, save all info to the HDF5 file
 *
 */

void DetectorWithAnalyticNonGaussianPSF::flushOutput()
{
    int Npixels    = 8;
    int Nsubpixels = 128;

    // Create the group in the HDF5 file.
    // We chose the same name as for DetectorWithMappedPSF

    hdf5File.createGroup("/PSF");

    // Generate and save the high resolution PSF (center of subfield)

    if (writeHighResolutionPSF)
    {
      Log.info("Writing high resolution PSF to the HDF5 file");

      arma::Mat<float> highResMap;
      makeHighResolutionPSF(highResMap, false, Npixels, Nsubpixels);
      hdf5File.writeArray("/PSF", "highResPSF", highResMap);
    }

    // Generate and save the diffused high resolution PSF (center of subfield)

    if (writeDiffusedPSF && includeChargeDiffusion)
    {
      Log.info("Writing diffused high resolution PSF to the HDF5 file");

      arma::Mat<float> highResDiffusedMap;
      makeHighResolutionPSF(highResDiffusedMap, true, Npixels, Nsubpixels);
      hdf5File.writeArray("/PSF", "diffusedPSF", highResDiffusedMap);
    }

    // Save the photometry info

    if (includePhotometry)
    {
        Log.info("Detector: Writing photometry to the HDF5 file");
	photometry.writePhotometry();
    }
    
}
