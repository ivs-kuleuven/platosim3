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

DetectorWithAnalyticNonGaussianPSF::DetectorWithAnalyticNonGaussianPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
: Detector(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure), sigma(nullptr)
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

    includePhotometry    = configParam.getBoolean("Photometry/IncludePhotometry");

    if (includePhotometry)
    {
        contaminationRadius = configParam.getInteger("Photometry/ContaminationRadius");                   // [pix]
        maskUpdateInterval  = configParam.getDouble("Photometry/MaskUpdateInterval") * 86400.;            // [s]                  
        filename            = configParam.getAbsoluteFilename("Photometry/TargetFileName");

        // Read and store the list of star IDs for which we want a lightcurve

        ifstream inputfile(filename);
        if (!inputfile) 
        {
            Log.error("DetectorWithAnalyticNonGaussianPSF::configure(): 'TargetFileName' file doesn't exist or is not readable: "  + filename);
            throw ConfigurationException("DetectorWithAnalyticNonGaussianPSF: wrong TargetFileName in configuration file");
        }

        photStarIDs.clear();
        while (getline(inputfile, line)) 
        {
            if (line == "" || line.find("#") == 0)
            {
                continue;
            }
            else
            {
                int starID = (unsigned int)stoi(line);
                photStarIDs.emplace_back(starID);
            }
        }

        // Initialize the maps containing the photometric information

        for (auto starID : photStarIDs)
        {
            inputFluxTarget[starID] = vector<double>(numExposures);  
            estimatedFluxTarget[starID] = vector<double>(numExposures); 
            varFluxTarget[starID] = vector<double>(numExposures);
            maskSizeTarget[starID] = vector<unsigned int>();  
            NSRtarget[starID] = vector<double>();   
            exposureNrOfMaskUpdate[starID] = vector<unsigned int>();
        }
    }


    // The configuration for the HDF5 contents
    
    writeFlatfieldMap = configParam.getBoolean("ControlHDF5Content/WriteFlatfieldMap");


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
    sigma->updateValue(time);
}











/**
 * \brief Interpolate and rotate PSF parameters and sum up all parts to calculate the intergal of the analytic PSF.
 * 
 * \param psf:        container to hold the result of the integration
 * \param x:          x position of the PSF
 * \param y:          y position of the PSF
 * \param r:          radial distance of the PSF to the optical axis
 * \param p:          azimuth angle of the PSF
 * \param Nsubpixels: number of subpixels per pixel (e.g. 128, set to 1 for no subpixels)
 **/

void DetectorWithAnalyticNonGaussianPSF::integrateAnalyticPSF(IntegralOfAnalyticSignalResponse& psf, double x, double y, double r, double p, int Nsubpixels)
{
    double ox = x - floor(x);
    double oy = y - floor(y);
    double s = (*sigma)() * Nsubpixels;
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
            double h  = (1. - w) * (*i)[c1 + 2] + w * (*i)[c2 + 2];
            double b  = s * ((1. - w) * (*i)[c1 + 3] + w * (*i)[c2 + 3]);
            double rr = s * ((1. - w) * (*i)[c1 + 4] + w * (*i)[c2 + 4]);
            double m  = (1. - w) * (*i)[c1 + 5] + w * (*i)[c2 + 5];
            double a  = (1. - w) * (*i)[c1 + 6] + w * (*i)[c2 + 6];

            psf.addPart(ox + pr * cos(p + pp), oy + pr * sin(p + pp), h, b, rr, m, p + a);
            psf.addPart(ox + pr * cos(p - pp), oy + pr * sin(p - pp), h, b, rr, m, p - a);
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

    Log.info("Detector: Integrating light for exposure " + to_string(exposureNr) + " with exposure time = " + to_string(exposureTime));

    integrateLight(exposureNr, startTime, exposureTime);

    // Include noise effects like readout noise, photon noise, full well saturation, etc.
    // Note: readOut() needs the exposure time to compute the open shutter smearing.

    Log.info("Detector: Adding noise effects to exposure " + to_string(exposureNr));

    readOut(exposureTime);

    // If photometric extraction was asked, apply it now

    if (includePhotometry)
    {
        Log.info("Detector: applying photometric extraction to exposure " + to_string(exposureNr));
        applyPhotometry(exposureNr);
    }

    // Write the CCD subfield, the bias map, and the smearing map to the HDF5 file

    Log.debug("Detector: Writing PixelMap, smearing map, and bias map #" + to_string(exposureNr) + " to HDF5 file.");

    writePixelMapsToHDF5(exposureNr);

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
    camera.exposeDetectorWithSkyBackground(*this, startTime, exposureTime, readoutTimeBeforeNextExposure);

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

    // Brighter-Fatter effect

    if(includeBFE)
    {
    		Log.debug("Detector: adding Brighter-Fatter effect");

       	applyBFE();
    }
    else
    {
        Log.debug("Detector: no Brighter-Fatter effect added");
    }
}







/**
 * \brief: Add the PSF of the star with given focal plane coordinates and flux level to the given map.
 *         As PSF we use an analytic non-Gaussian function.
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

    double ccdOrientation = getOrientationAngle();
    p -= ccdOrientation;

    integrateAnalyticPSF(psf, column0, row0, r, p);

    for (int y = max(0, sy); y < min((int)map.n_rows, sy + size); y++)
        for (int x = max(0, sx); x < min((int)map.n_cols, sx + size); x++)
            map.at(y, x) += psf(x - sx, y - sy) * flux;

    return true;
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
 * \brief: Create a high-resolution map of the PSF in the center of the subfield
 *         
 * \param map         The PSF will be written in this 2D array, the size need not be allocated.
 * \param Npixels     Size of the high-res map in pixels
 * \param Nsubpixels  (sqrt of the) number of subpixels per pixel 
 *
 * \return           'map' will be modified to contain the high-resolution PSF
 *                   Its size will be  (Npixels * Nsubpixels) x  (Npixels * Nsubpixels)
 */

void DetectorWithAnalyticNonGaussianPSF::makeHighResolutionPSF(arma::Mat<float> &highResMap, int Npixels, int Nsubpixels)
{
    // Put the PSF right in the middle of the (high-res) (sub)pixel map

    double row0 = Npixels / 2.0;
    double column0 = Npixels / 2.0;

    // The high-res (sub)pixel map will be placed in the middle of the regular subfield.
    // Derive the focal plane coordinates of the middle of the subfield. These are needed
    // to later on derive the angular distance from the optical axis.

    double middleRowSubfield = subFieldZeroPointRow + numRowsPixelMap / 2.0;
    double middleColSubfield = subFieldZeroPointColumn +  numColumnsPixelMap / 2.0;
    double xFP, yFP;
    tie(xFP, yFP) = pixelToFocalPlaneCoordinates(middleRowSubfield, middleColSubfield);


    double s = (*sigma)();
    double diffusionKernelWidth = 0.;

    if (includeChargeDiffusion) 
    {
        diffusionKernelWidth = chargeDiffusionStrength;
        s = sqrt(s * s + diffusionKernelWidth * diffusionKernelWidth);
    }

    int size = Npixels * Nsubpixels;
    highResMap.set_size(size, size);

    int sx = (int)floor(column0*Nsubpixels - (size - 1.) / 2.);
    int sy = (int)floor(row0*Nsubpixels - (size - 1.) / 2.);

    // Construct the PSF around the central pixel coordinates

    IntegralOfAnalyticSignalResponse psf(size, diffusionKernelWidth*Nsubpixels);
    double r = rad2deg(camera.getGnomonicRadialDistanceFromOpticalAxis(xFP, yFP));
    double p = atan2(yFP, xFP);

    double ccdOrientation = getOrientationAngle();
    p -= ccdOrientation;

    integrateAnalyticPSF(psf, column0, row0, r, p, Nsubpixels);

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







/**
 *  \brief Before destroying this object, save all info to the HDF5 file
 * 
 */ 

void DetectorWithAnalyticNonGaussianPSF::flushOutput()
{
    int Npixels = 8;
    int Nsubpixels = 128;

    // Create the group in the HDF5 file. We chose the same name as for DetectorWithMappedPSF

    hdf5File.createGroup("/PSF");
    
    // Generate the high-resolution map

    arma::Mat<float> highResMap;
    makeHighResolutionPSF(highResMap, Npixels, Nsubpixels);

    // Save the map to HDF5

    hdf5File.writeArray("/PSF", "HighResPSFmapCenterSubfield", highResMap);

    // Save the photometry info

    if (includePhotometry)
    {
        Log.info("Writing photometry to the HDF5 file");

        hdf5File.createGroup("/Photometry");
        hdf5File.createGroup("/Photometry/Masks");
        hdf5File.createGroup("/Photometry/Lightcurves");

        string groupName = "/Photometry/Masks";
        string arrayName = "exposureNrOfMaskUpdate";
        int starID = photStarIDs[0];    // Doesn't matter which one we take the update exposure Nrs are the same for all targets.
        hdf5File.writeArray(groupName, arrayName, exposureNrOfMaskUpdate[starID].data(), exposureNrOfMaskUpdate[starID].size());

        for (auto starID : photStarIDs)
        {
            string starName = to_string(starID);
            groupName = "/Photometry/Lightcurves/starID" + starName;
            hdf5File.createGroup(groupName);

            arrayName = "inputFlux";
            hdf5File.writeArray(groupName, arrayName, inputFluxTarget[starID].data(), inputFluxTarget[starID].size());

            arrayName = "estimatedFlux";
            hdf5File.writeArray(groupName, arrayName, estimatedFluxTarget[starID].data(), estimatedFluxTarget[starID].size());

            groupName = "/Photometry/Masks/starID" + starName;
            hdf5File.createGroup(groupName);

            arrayName = "maskSize";
            hdf5File.writeArray(groupName, arrayName, maskSizeTarget[starID].data(), maskSizeTarget[starID].size());

            arrayName = "maskNSR";
            hdf5File.writeArray(groupName, arrayName, NSRtarget[starID].data(), NSRtarget[starID].size());

            for(auto iter = rowIndexOfMaskOfTarget[starID].begin(); iter != rowIndexOfMaskOfTarget[starID].end(); ++iter)
            {
                const unsigned int exposureNumber = iter->first;
                
                stringstream myStream;
                myStream << "Exposure" << setfill('0') << setw(6) << exposureNumber;
                groupName = "/Photometry/Masks/starID" + starName + "/" + myStream.str();
                hdf5File.createGroup(groupName);

                arrayName = "maskRowIndices"; 
                hdf5File.writeArray(groupName, arrayName, rowIndexOfMaskOfTarget[starID][exposureNumber].data(), rowIndexOfMaskOfTarget[starID][exposureNumber].size());
            
                arrayName = "maskColumnIndices";
                hdf5File.writeArray(groupName, arrayName, colIndexOfMaskOfTarget[starID][exposureNumber].data(), colIndexOfMaskOfTarget[starID][exposureNumber].size());
            }
        }
    } // end if includePhotometry
} // end flushOutput()











/**
 * \brief Extract the photometric light curve for a specified list of stars
 *
 * TODO: - better error catching when the stars for which a lightcurve is requested are (sometimes) not in the subfield
 *       - better treatment when there are no contaminants
 */

void DetectorWithAnalyticNonGaussianPSF::applyPhotometry(const unsigned int exposureNr)
{
    const unsigned int zeroBasedExposureNr = exposureNr - beginExposureNr;

    const double varianceRON = sqrt(pow(readoutNoise, 2) + pow(frontEndElectronics->getReadoutNoise(), 2));      // [electrons / pixel]
 
    // Make a (deep) copy of the pixelMap on which we can do some reductions without altering the original pixelMap

    arma::Mat<float> image(pixelMap);

    // Subtract the bias

    float meanBias = 0.0;
    if (subFieldZeroPointColumn <  numColumns / 2)
    {
        meanBias = arma::mean(arma::mean(biasMapLeft));
    }
    else
    {
        meanBias = arma::mean(arma::mean(biasMapRight));
    }
    image -= meanBias;

    // Correct for open-shutter smearing

    image.each_row() -= arma::mean(smearingMap - meanBias, 0);

    // Convert from [ADU] to [electrons] using the gain

    if (subFieldZeroPointColumn <  numColumns / 2)
    {
        image /= combinedGainLeft;
    }
    else
    {
        image /= combinedGainRight;
    }
    
    // Subtract the sky background

    const double skyBackground = camera. getTotalSkyBackground();                // [photons/pixel/exposure]
    image -= throughputMap * skyBackground;                                      // [e-/pixel/exposure]


    // Loop over all targets for which you need a lightcurve

    const int Ntargets = photStarIDs.size();                                     // Nr of stars for which we want a lightcurve

    for (int n = 0; n < Ntargets; n++)
    {
        // Collect info on the position and the input flux of the target

        int starID = photStarIDs[n];

        double time;                                                      // Time stamp of the last exposure         [s]
        double xFPtarget;                                                 // Mean x-coordinate in the focal plane    [mm]
        double yFPtarget;                                                 // Mean y-coordinate in the focal plane    [mm] 
        double rowTarget;                                                 // Mean row coordinate in the subfield     [pix] 
        double colTarget;                                                 // Mean column coordinate in the subfield  [pix]
        double fluxTarget;                                                // Total flux during the exposure          [photons/exposure]
    
        tie(time, xFPtarget, yFPtarget, rowTarget, colTarget, fluxTarget) = camera.getInfoForTheMostRecentExposureForStar(starID);

        inputFluxTarget.at(starID).at(zeroBasedExposureNr) = fluxTarget;

        // If this is the first exposure, or it's already 2 weeks ago that the mask was updated,
        // generate the mask of the current target. 

        double timeSinceLastMaskUpdate = 0.0;

        if (exposureNr != beginExposureNr)
        {
            timeSinceLastMaskUpdate = (exposureNr - exposureNrOfMaskUpdate.at(starID).back()) * cycleTime;  // [s]
        }

        if ((exposureNr == beginExposureNr) || (timeSinceLastMaskUpdate > maskUpdateInterval))
        {
            Log.debug("Detector::applyPhotometry: updating mask of star ID " + to_string(starID) + " for exposure " + to_string(exposureNr));
            Log.debug("Detector::applyPhotometry: creating single-target and contamination maps");

            arma::Mat<float> singleTargetMap(numRowsPixelMap, numColumnsPixelMap);
            arma::Mat<float> contaminantMap(numRowsPixelMap, numColumnsPixelMap);

            // Create a noiseless subfield as if there was only the flux of this single target

            singleTargetMap.zeros();
            double r = rad2deg(camera.getGnomonicRadialDistanceFromOpticalAxis(xFPtarget, yFPtarget));
            double p = atan2(yFPtarget, xFPtarget);
            bool success = addFluxToMap(singleTargetMap, rowTarget, colTarget, r, p, fluxTarget);

            // Create a noiseless subfield of only the possible contaminants

            contaminantMap.zeros();

            int numContaminants = 0;
            starInfoIterator begin, end;
            tie(begin, end) = camera.getInfoForTheMostRecentExposureForAllStars();
            for (auto it = begin; it != end; it++)
            {
                if (it->first == starID) continue;                        // A star is never its own contaminant 
                double xFPcont =  (it->second)[0] / (it->second)[5];      // [mm]
                double yFPcont =  (it->second)[1] / (it->second)[5];      // [mm]
                double rowCont =  (it->second)[2] / (it->second)[5];      // [pix]
                double colCont =  (it->second)[3] / (it->second)[5];      // [pix]
                double fluxCont = (it->second)[4];                        // [photons/exposure]

                Log.debug(to_string(it->first) + ": " + to_string(rowCont) + ", " + to_string(colCont) + ", " + to_string(fluxCont));

                // Skip the contaminants that are too distant from the target to have any effect

                if ((abs(colCont - colTarget) > contaminationRadius) or (abs(rowCont - rowTarget) > contaminationRadius))
                    continue;
                
                // Add the PSF of the contaminant to the contaminant map
                
                r = rad2deg(camera.getGnomonicRadialDistanceFromOpticalAxis(xFPcont, yFPcont));
                p = atan2(yFPcont, xFPcont);
                success = addFluxToMap(contaminantMap, rowCont, colCont, r, p, fluxCont);
                numContaminants++;
            }

            Log.debug("Detector::applyPhotometry: Found " + to_string(numContaminants) + " contaminants for star ID " + to_string(starID));
            Log.debug("Detector::applyPhotometry: selecting which pixels belong to the mask");

            // For the mask of our target will only consider a 7x7 area around the barycenter. Get the boundaries of that area.

            const int minRow = max(0, int(rowTarget)-3);
            const int maxRow = min(int(numRowsPixelMap) - 1, int(rowTarget)+3);                         // maxRow inclusive
            const int minCol = max(0, int(colTarget)-3);
            const int maxCol = min(int(numColumnsPixelMap) - 1, int(colTarget)+3);                      // maxCol inclusive
                
            // For the pixels in the designated area around our target, compute the variance and the noise/signal ratio of the signal.

            arma::Mat<float> NSRmap(numRowsPixelMap, numColumnsPixelMap, arma::fill::zeros); 
            arma::Mat<float> varianceMap(numRowsPixelMap, numColumnsPixelMap, arma::fill::zeros);

            vector<double> flatNSRmap;   
            for (int irow = minRow; irow <= maxRow; irow++)
            {
                for (int icol = minCol; icol <= maxCol; icol++)
                {
                    // We assume photon noise, so the variance equals the flux. We multiply by the throughput so that both terms
                    // are expressed in [e-/exposure]. 

                    varianceMap(irow, icol) = (singleTargetMap(irow, icol) + contaminantMap(irow, icol) + skyBackground) * throughputMap(irow, icol) + varianceRON ;
                    NSRmap(irow, icol) = sqrt(varianceMap(irow, icol)) / singleTargetMap(irow, icol); 
                    flatNSRmap.push_back(NSRmap(irow, icol));
                }
            }

            // Order the pixels in the (flattened) NSR map from low to high N/S ratio (i.e. high to low S/N)

            vector<unsigned int> indices(flatNSRmap.size());
            iota(indices.begin(), indices.end(), 0);
            stable_sort(indices.begin(), indices.end(), [&flatNSRmap](unsigned int i, unsigned int j) {return flatNSRmap[i] < flatNSRmap[j];});
            vector<unsigned int> rowIndex(flatNSRmap.size());           
            vector<unsigned int> colIndex(flatNSRmap.size());
            const int N = maxRow-minRow +1; 
            for (int i = 0; i < rowIndex.size(); i++)          // Transform from indices in flatNSRmap to indices in NSRmap
            {
                rowIndex[i] = minRow + (unsigned int)(indices[i]) / N;
                colIndex[i] = minCol + (unsigned int)(indices[i]) % N; 
            }

            // Build the mask, starting with the pixel with the best NSR, adding one pixel at the time,
            // with the condition that adding a pixel should contribute more to the aggregated signal than to the aggregated noise.

            // Initialize with the first pixel

            double aggregatedVariance            = varianceMap(rowIndex[0], colIndex[0]);
            double aggregatedSingleTargetFlux    = singleTargetMap(rowIndex[0], colIndex[0]);
            double aggregatedObservedTargetFlux  = image(rowIndex[0], colIndex[0]);
            double aggregatedNSR                 = NSRmap(rowIndex[0], colIndex[0]);
            maskSizeTarget[starID].push_back(1);
 
            rowIndexOfMaskOfTarget[starID][exposureNr] = {rowIndex[0]}; 
            colIndexOfMaskOfTarget[starID][exposureNr] = {colIndex[0]};

            // Then add other pixels

            for (int i = 1; i < rowIndex.size(); i++)
            {
                double temp = sqrt(aggregatedVariance + varianceMap(rowIndex[i], colIndex[i])) / (aggregatedSingleTargetFlux + singleTargetMap(rowIndex[i], colIndex[i]));
                if (temp < aggregatedNSR)
                {
                    // The aggregated Noise / Signal ratio improved by adding a pixel, so include the pixel in the mask

                    aggregatedVariance           += varianceMap(rowIndex[i], colIndex[i]);
                    aggregatedSingleTargetFlux   += singleTargetMap(rowIndex[i], colIndex[i]);
                    aggregatedObservedTargetFlux += image(rowIndex[i], colIndex[i]);
                    aggregatedNSR = temp;
                    maskSizeTarget.at(starID).back() += 1;
                    rowIndexOfMaskOfTarget.at(starID).at(exposureNr).push_back(rowIndex[i]);
                    colIndexOfMaskOfTarget.at(starID).at(exposureNr).push_back(colIndex[i]);
                }
                else
                {
                    // The aggregated Noise/Signal ratio did not improve by adding this pixel. Not only can we ignore exclude this pixel from the 
                    // mask, but also all subsequent ones that have an even worse noise/signal ratio. So finalize the mask for this target, and 
                    // then break out of the for-loop.

                    estimatedFluxTarget.at(starID).at(zeroBasedExposureNr) = aggregatedObservedTargetFlux; 
                    varFluxTarget.at(starID).at(zeroBasedExposureNr) = aggregatedVariance; 
                    NSRtarget.at(starID).push_back(aggregatedNSR);

                    // Disregard all other pixels of the window around the target star: they all contribute more to the noise than to the signal.

                    break;
                }
            }

            // Update the exposure nr of this mask update to the current exposure number

            exposureNrOfMaskUpdate.at(starID).push_back(exposureNr);

        }
        else
        {
            // For all other exposures, simply use the same (most recent) mask. We reuse the NSR of the previous mask, so no need to recompute it.

            Log.debug("Detector::applyPhotometry: extracting flux for target ID " + to_string(starID) + " for exposure " + to_string(exposureNr) + " with an old mask");

            const unsigned int exposureNrOfLastMaskUpdate =  exposureNrOfMaskUpdate.at(starID).back();
            estimatedFluxTarget.at(starID).at(zeroBasedExposureNr) = 0.0;
            for (int j = 0; j < maskSizeTarget[starID].back(); j++)
            {
                const unsigned int rowIndex = rowIndexOfMaskOfTarget.at(starID).at(exposureNrOfLastMaskUpdate).at(j);
                const unsigned int colIndex = colIndexOfMaskOfTarget.at(starID).at(exposureNrOfLastMaskUpdate).at(j);
                estimatedFluxTarget.at(starID).at(zeroBasedExposureNr) += image(rowIndex, colIndex);
            }
        }
    } // end loop over all targets for which we want light curves
} // end applyPhotometry()
