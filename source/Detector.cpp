#include "Detector.h"

/**
 * \brief Add a part to the definite integral of the analytic PSF and adjust the normalization factor accordingly.
 *
 * \param ox:    relative x offset
 * \param oy:    relative y offset
 * \param h:     relative height
 * \param sigma: width of Gaussian term
 * \param r:     strength and type of periodic term
 * \param rho:   distance between Gaussian term and periodic term
 * \param phi:   orientation of periodic term
 *
 * \return reference to itself
 **/

IntegralOfAnalyticSignalResponse& IntegralOfAnalyticSignalResponse::addPart(double ox, double oy, double h, double sigma, double r, double rho, double phi)
{
    using Faddeeva::erf;

    ox += (size - (size & 1)) / 2.;
    oy += (size - (size & 1)) / 2.;

    erfxr.emplace_back(size + 1);
    erfyr.emplace_back(size + 1);
    double sr = 1. / sqrt(2.) / sigma;
    double fr1 = sqrt(M_PI * fabs(h) * sigma * sigma / (r != 0.? 4.: 2.));
    double fr2 = h < 0.? -fr1: fr1;

    // dsigma is the Gaussian diffusion kernel width

    if (dsigma != 0.)
        sr /= sqrt(2. * sr * sr * dsigma * dsigma + 1.);

    erfxr.back()[0] = erf(-sr * ox);
    erfyr.back()[0] = erf(-sr * oy);
    for (unsigned i = 0; i < size; i++)
    {
        erfxr.back()[i] = ((erfxr.back()[i + 1] = erf(sr * (i + 1. - ox))) - erfxr.back()[i]) * fr1;
        erfyr.back()[i] = ((erfyr.back()[i + 1] = erf(sr * (i + 1. - oy))) - erfyr.back()[i]) * fr2;
    }
    n += 4. * fr1 * fr2;

    if (r != 0.)
    {
        erfxc.emplace_back(size + 1);
        erfyc.emplace_back(size + 1);
        double delta = 2. * M_PI * sigma * sigma / r / r;
        complex<double> sc = 1. / sqrt(2.) / sigma * sqrt(complex<double>(1., delta));
        complex<double> xc = complex<double>(0., M_PI / fabs(r) * sqrt(rho) * cos(phi)) / sc;
        complex<double> yc = complex<double>(0., M_PI / fabs(r) * sqrt(rho) * sin(phi)) / sc;
        complex<double> fc = sqrt((r < 0.? -fr1: fr1) * fr2 * exp(-M_PI * rho / (1. + delta * delta) * complex<double>(delta, 1.)) / complex<double>(1., delta));

        if (dsigma != 0.) {
            complex<double> dc = sqrt(2. * sc * sc * dsigma * dsigma + 1.);
            sc /= dc;
            xc /= dc;
            yc /= dc;
        }

        erfxc.back()[0] = erf(xc - sc * ox);
        erfyc.back()[0] = erf(yc - sc * oy);
        for (unsigned i = 0; i < size; i++)
        {
            erfxc.back()[i] = ((erfxc.back()[i + 1] = erf(xc + sc * (i + 1. - ox))) - erfxc.back()[i]) * fc;
            erfyc.back()[i] = ((erfyc.back()[i + 1] = erf(yc + sc * (i + 1. - oy))) - erfyc.back()[i]) * fc;
        }
        n -= 4. * (fc.real() * fc.real() - fc.imag() * fc.imag());
    }
    return *this;
}









/**
 * \brief Get integral of analytic PSF for a (sub)pixel.
 *
 * \param i:    x direction
 * \param j:    y direction
 * \param norm: apply normalization
 *
 * \return normalized or unnormalized integral of analytic PSF for (sub)pixel (i, j)
 **/

double IntegralOfAnalyticSignalResponse::operator()(unsigned i, unsigned j, bool norm)
{
    if (i >= size || j >= size)
    {
        return 0.;
    }

    double ret = 0.;
    for (unsigned k = 0; k < erfxr.size() && k < erfyr.size(); k++)
    {
        ret += erfxr[k][i] * erfyr[k][j];
    }

    for (unsigned k = 0; k < erfxc.size() && k < erfyc.size(); k++)
    {
        ret -= erfxc[k][i].real() * erfyc[k][j].real() - erfxc[k][i].imag() * erfyc[k][j].imag();
    }

    return norm? ret / n: ret;
}








/**
 * \brief Constructor.
 * 
 * \details
 * 
 * The constructor initializes the groups in the HDF5 file where the different maps (i.e. pixel map,
 * bias register map, smearing map, etc.) will be saved. 
 * 
 * The following maps are initialized:
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

Detector::Detector(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
: HDF5Writer(hdf5file),
  includeCosmicsInSubField(true), includeCosmicsInSmearingMap(true), includeCosmicsInBiasMap(true),
  includeBFE(true),
  includeDarkSignal(true),
  includePhotonNoise(true),
  includeReadoutNoise(true),
  includeCTIeffects(true), 
  includeOpenShutterSmearing(true), 
  includeQuantumEfficiency(true),
  includeNaturalVignetting(true),
  includeMechanicalVignetting(true),
  includeParticulateContamination(true),
  includeMolecularContamination(true),
  includeFullWellSaturation(true),
  includeDigitalSaturation(true),
  internalTime(0.0), camera(camera),
  temperatureGenerator(detectorTemperatureGenerator)
{
    // Parse the parameters from the configuration file.

    configure(configParam);

    this->readoutTimeBeforeNextExposure = readoutTimeBeforeNextExposure;
    this->readoutTimeDuringNextExposure = readoutTimeDuringNextExposure;

    // Create the groups in the HDF5 file where the different maps (i.e. pixel map,
    // bias register map, smearing map, etc.) will be saved. This needs to be done
    // BEFORE other methods write arrays to HDF5.

    initHDF5Groups();

    // Front-end electronics

    frontEndElectronics = new FrontEndElectronics(configParam, hdf5file, feeTemperatureGenerator);

    // Allocate memory for the different maps

    pixelMap.zeros(numRowsPixelMap, numColumnsPixelMap);
    biasMapLeft.zeros(numRowsBiasMap, numColumnsBiasMap);
    biasMapRight.zeros(numRowsBiasMap, numColumnsBiasMap);

    smearingMap.zeros(numRowsSmearingMap, numColumnsPixelMap);
    throughputMap.ones(numRowsPixelMap, numColumnsPixelMap);

    // If we are going to apply open-shutter smearing, we have to know which pixels are within
    // the FOV (relevant only in case of mechanical vignetting).  When mechanical vignetting is
    // disabled, all pixels of the detector are inside the FOV.

    if(includeOpenShutterSmearing)
    {
        // Mechanical vignetting map:
        //  - no mechanical vignetting: all pixels of the sub-field inside FOV -> all values set to one
        //  - mechanical vignetting: set value of the pixels in the sub-field outside FOV to zero (others should be one) -> on creation of the throughput map

        mechanicalVignettingMask.ones(numRowsPixelMap, numColumnsPixelMap);

        // Number of exposed rows in each column:
        // - no mechanical vignetting: all exposed rows inside FOV (numRows - firstRowExposed)
        // - mechanical vignetting: count the exposed rows (i.e. from firstRowExposed) that are inside FOV

        numExposedRowsInFOV.zeros(numColumnsPixelMap);

        if(!includeMechanicalVignetting)
            numExposedRowsInFOV.fill(numRows - firstRowExposed);
    }

    // Check whether the gain values for left- and right-hand side of the CCD are not too far apart

    checkGain();

    // Set the seeds of the random number generators

    darkSignalGenerator.seed(darkSignalSeed);
    darkNoiseGenerator.seed(darkSignalSeed + 1);

    photonNoiseGenerator.seed(photonNoiseSeed);
    readoutNoiseGenerator.seed(readoutNoiseSeed);

    cosmicHitRateGenerator.seed(cosmicSeed);
    cosmicEntryRowGenerator.seed(cosmicSeed + 1);
    cosmicEntryColumnGenerator.seed(cosmicSeed + 2);
    cosmicEntryAngleGenerator.seed(cosmicSeed + 3);
    cosmicTrailLengthGenerator.seed(cosmicSeed + 4);
    cosmicIntensityGenerator.seed(cosmicSeed + 5);
    decimalNumCosmicHitsGenerator.seed(cosmicSeed + 6);

    decimalNumCosmicHitsDistribution = uniform_real_distribution<double>(0, 1);

    firstExposure = true;
}








/**
 * Destructor.
 *
 */
Detector::~Detector()
{
    flushOutput();

    delete frontEndElectronics;
}







/**
 * \brief Update the time dependent parameters of the Detector to their 
 *        value at the given time point
 *
 * \param time: current time
 */

void Detector::updateParameters(double time)
{

}











/**
 * \brief Configure the Detector object using the ConfigurationParameters
 * 
 * \param configParam: Configuration parameters
 */

 void Detector::configure(ConfigurationParameters &configParam)
 {
    // Configuration parameters for the CCD detector

    string ccdPosition                  = configParam.getString("CCD/Position");

    if (ccdPosition == "Custom")
    {
        originOffsetX         = configParam.getDouble("CCD/OriginOffsetX");        // [mm]
        originOffsetY         = configParam.getDouble("CCD/OriginOffsetY");        // [mm]
        orientationAngle      = deg2rad(configParam.getDouble("CCD/Orientation")); // [rad]
        numRows               = configParam.getInteger("CCD/NumRows");             // [pixels]
        numColumns            = configParam.getInteger("CCD/NumColumns");          // [pixels]
        firstRowExposed       = configParam.getInteger("CCD/FirstRowExposed");     // [pixels]
    }
    else
    {
        int idx = stoi(ccdPosition) - 1;  // Positions are named [1, 2, 3, 4] while the index into vector starts at 0

        originOffsetX         = configParam.getDoubleAt("CCDPositions/OriginOffsetX", idx);
        originOffsetY         = configParam.getDoubleAt("CCDPositions/OriginOffsetY", idx);
        orientationAngle      = deg2rad(configParam.getDoubleAt("CCDPositions/Orientation", idx));
        numRows               = configParam.getIntegerAt("CCDPositions/NumRows", idx);
        numColumns            = configParam.getIntegerAt("CCDPositions/NumColumns", idx);

        isFastCamera          = configParam.getString("Telescope/GroupID") == "Fast";
        
        if (isFastCamera)
        {
            firstRowExposed       = configParam.getIntegerAt("CCDPositions/FirstRowForFastCamera", idx);
        }
        else
        {
            firstRowExposed       = configParam.getIntegerAt("CCDPositions/FirstRowForNormalCamera", idx);
        }
    }

    Log.debug("Detector: selected ccdPosition = " + ccdPosition);
    Log.debug("Detector: CCD originOffsetX, originOffsetY = " + to_string(originOffsetX) + ", " + to_string(originOffsetY) + " mm");
    Log.debug("Detector: CCD orientationAngle = " + to_string(rad2deg(orientationAngle)) + " deg");
    Log.debug("Detector: CCD numRows, numColumns, firstRow = " + to_string(numRows) + ", " + to_string(numColumns) + ", " + to_string(firstRowExposed));

    pixelSize                           = configParam.getDouble("CCD/PixelSize");
//    quantumEfficiency                   = configParam.getDouble("CCD/QuantumEfficiency/Efficiency");                  // FIXME: No commented out lines of code. To be removed or not?
//    refAngleQE                          = configParam.getDouble("CCD/QuantumEfficiency/RefAngle");
//    relativeRefEfficiencyQE             = configParam.getDouble("CCD/QuantumEfficiency/RelativeRefEfficiency");
    meanQE                              = configParam.getDouble("CCD/QuantumEfficiency/MeanQuantumEfficiency");
    meanAngleDependencyQE               = configParam.getDouble("CCD/QuantumEfficiency/MeanAngleDependency");
//    expectedValueQuantumEfficiency      = configParam.getDouble("CCD/QuantumEfficiency/ExpectedValue");
    includeCosmicsInSubField            = configParam.getBoolean("Sky/IncludeCosmicsInSubField");
    includeCosmicsInSmearingMap         = configParam.getBoolean("Sky/IncludeCosmicsInSmearingMap");
    includeCosmicsInBiasMap             = configParam.getBoolean("Sky/IncludeCosmicsInBiasMap");
    cosmicHitRate                       = configParam.getDouble("Sky/Cosmics/CosmicHitRate");
    cosmicTrailLength                   = configParam.getDoubleVector("Sky/Cosmics/TrailLength");
    cosmicIntensity                     = configParam.getDoubleVector("Sky/Cosmics/Intensity");
    darkCurrent                         = configParam.getDouble("CCD/DarkSignal/DarkCurrent");
    dsnu                                = configParam.getDouble("CCD/DarkSignal/DSNU");
    darkCurrentStability                = configParam.getDouble("CCD/DarkSignal/Stability");
    includeBFE                          = configParam.getBoolean("CCD/IncludeBFE");
    rangeBFE                            = configParam.getInteger("CCD/BFE/Range");
    p0BFE                               = configParam.getDouble("CCD/BFE/p0");
    p1BFE                               = configParam.getDouble("CCD/BFE/p1");
    refFluxBFE                          = configParam.getDouble("CCD/BFE/RefFlux");

    fullWellSaturationLimit             = configParam.getLong("CCD/FullWellSaturation");
    digitalSaturationLimit              = configParam.getLong("CCD/DigitalSaturation");
    readoutNoise                        = configParam.getDouble("CCD/ReadoutNoise");
    expectedValueNaturalVignetting      = configParam.getDouble("CCD/Vignetting/NaturalVignetting/ExpectedValue");
    radiusFOV                           = deg2rad(configParam.getDouble("CCD/Vignetting/MechanicalVignetting/RadiusFOV"));
    particulateContaminationEfficiency  = configParam.getDouble("CCD/Contamination/ParticulateContaminationEfficiency");
    molecularContaminationEfficiency    = configParam.getDouble("CCD/Contamination/MolecularContaminationEfficiency");

    refValueGainLeft          = configParam.getDouble("CCD/Gain/RefValueLeft");
    refValueGainRight         = configParam.getDouble("CCD/Gain/RefValueRight");
    gainStability             = configParam.getDouble("CCD/Gain/Stability");
    gainAllowedDifference     = configParam.getDouble("CCD/Gain/AllowedDifference");

    readoutMode               = configParam.getString("CCD/ReadoutMode/ReadoutMode");

    if(readoutMode == "Partial")
    {
    	firstRowPartialReadout = configParam.getInteger("CCD/ReadoutMode/Partial/FirstRowReadout");
    	numRowsPartialReadout = configParam.getInteger("CCD/ReadoutMode/Partial/NumRowsReadout");
    }
    else if (readoutMode != "Nominal")
    {
    	Log.error("Detector::configure(): Unknown readout mode specification in configuration file: "  + readoutMode);
    	throw ConfigurationException("Detector: Unknown readout mode specification in configuration file");
    }

    serialTransferTime = configParam.getDouble("CCD/SerialTransferTime") * 1E-9;			  // [ns] -> [s]
    parallelTransferTime = configParam.getDouble("CCD/ParallelTransferTime") * 1E-6;		  // [µs] -> [s]
    parallelTransferTimeFast = configParam.getDouble("CCD/ParallelTransferTimeFast") * 1E-6;  // [µs] -> [s]

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

//    polarizationEfficiency          = configParam.getDouble("CCD/Polarization/Efficiency");
//    refAnglePolarization            = configParam.getDouble("CCD/Polarization/RefAngle");
    expectedValuePolarization       = configParam.getDouble("CCD/Polarization/ExpectedValue");

    nominalOperatingTemperature     = configParam.getDouble("CCD/NominalOperatingTemperature");

    includeParticulateContamination = configParam.getBoolean("CCD/IncludeParticulateContamination");
    includeMolecularContamination   = configParam.getBoolean("CCD/IncludeMolecularContamination");
    includeDarkSignal               = configParam.getBoolean("CCD/IncludeDarkSignal");
    includePhotonNoise              = configParam.getBoolean("CCD/IncludePhotonNoise");
    includeReadoutNoise             = configParam.getBoolean("CCD/IncludeReadoutNoise");
    includeCTIeffects               = configParam.getBoolean("CCD/IncludeCTIeffects");
    includeOpenShutterSmearing      = configParam.getBoolean("CCD/IncludeOpenShutterSmearing");
    includeQuantumEfficiency        = configParam.getBoolean("CCD/IncludeQuantumEfficiency");
    includeNaturalVignetting        = configParam.getBoolean("CCD/IncludeNaturalVignetting");
    includeMechanicalVignetting     = configParam.getBoolean("CCD/IncludeMechanicalVignetting");
    includePolarization             = configParam.getBoolean("CCD/IncludePolarization");
    includeFullWellSaturation       = configParam.getBoolean("CCD/IncludeFullWellSaturation");
    includeDigitalSaturation        = configParam.getBoolean("CCD/IncludeDigitalSaturation");
    includeQuantisation             = configParam.getBoolean("CCD/IncludeQuantisation");

    // Configuration parameters for the subfield

    subFieldZeroPointRow    = configParam.getInteger("SubField/ZeroPointRow");
    subFieldZeroPointColumn = configParam.getInteger("SubField/ZeroPointColumn");
    numRowsPixelMap         = configParam.getInteger("SubField/NumRows");
    numColumnsPixelMap      = configParam.getInteger("SubField/NumColumns");
    numRowsBiasMap          = configParam.getInteger("SubField/NumBiasPrescanRows");
    numColumnsBiasMap       = configParam.getInteger("SubField/NumBiasPrescanColumns");
    numRowsSmearingMap      = configParam.getInteger("SubField/NumSmearingOverscanRows");

    Log.debug("Detector: Subfield zero point (row, col) = (" + to_string(subFieldZeroPointRow) + ", " + to_string(subFieldZeroPointColumn) + ")");
    Log.debug("Detector: Subfield center point (row, col) = (" + to_string(subFieldZeroPointRow + numRowsPixelMap/2) 
                                                               + ", " + to_string(subFieldZeroPointColumn + numColumnsPixelMap/2) + ")");
    Log.debug("Detector: Subfield nr of rows = " + to_string(numRowsPixelMap));
    Log.debug("Detector: Subfield nr of columns = " + to_string(numColumnsPixelMap));

    // No parallel over-scan in case of partial readout

    if(readoutMode == "Partial")
    {
    	Log.info("No smearing map for partial readout");
    	numRowsSmearingMap = 0;
    }

    // Configuration parameters for the noise source random seeds

    readoutNoiseSeed        = configParam.getLong("RandomSeeds/ReadOutNoiseSeed");
    photonNoiseSeed         = configParam.getLong("RandomSeeds/PhotonNoiseSeed");
    cosmicSeed              = configParam.getLong("RandomSeeds/CosmicSeed");
    darkSignalSeed          = configParam.getLong("RandomSeeds/DarkSignalSeed");


    // Get the sequential number of the very first exposure

    beginExposureNr         = configParam.getInteger("ObservingParameters/BeginExposureNr");

    numEdgePixels = 0;

    sendImagettesToClient = false;
 }










/**
 * \brief: Take an exposure with the detector starting at the given time.
 *         The light is integrated during the given exposure time, during which 
 *         the detector experiences the effects of jitter and thermo-elastic telescope 
 *         drift. The background is assumed uniform for the whole subfield.
 *         Afterwards, the collected light is read out, convolving the image with the
 *         point spread function and adding various noise effects.
 *
 * \param exposureNr:   sequential number of the exposure
 * \param startTime:    Starting time of the exposure [s].
 * \param exposureTime: Duration of the exposure [s].
 * 
 * \return endTime:     Time after the exposure (startTime + exposureTime + readoutTime)
 *
 * \pre Sub-pixel, pixel, bias register, and smearing map filled with values from previous exposure.
 *
 * \post Pixel unit in the pixel, bias register, and smearing maps: [ADU]
 */

double Detector::takeExposure(int exposureNr, double startTime, double exposureTime)
{
    // Advance the internal clock until the given start time

    internalTime = startTime;

    // check if there are new updates for the window position

    if (getWinPositionFromServer)
    {
        setWinPosition();
    }

    // Integration of point sources and background, taking into account jitter + drift.

    Log.info("Detector: Integrating light for exposure " + to_string(exposureNr) + " with exposure time = " + to_string(exposureTime));

    integrateLight(exposureNr, startTime, exposureTime);

    // Include noise effects like readout noise, photon noise, full well saturation, etc.
    // Note: readOut() needs the exposure time to compute the open shutter smearing.

    Log.info("Detector: Adding noise effects to exposure " + to_string(exposureNr));

    readOut(exposureTime);

    // Write the CCD subfield, the bias map, and the smearing map to the HDF5 file

    Log.debug("Detector: Writing PixelMap, smearing map, bias map and throughputMap #" + to_string(exposureNr) + " to HDF5 file.");

    writePixelMapsToHDF5(exposureNr);

    // Advance the internal clock

    internalTime += exposureTime + readoutTimeBeforeNextExposure;

    return internalTime;
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

    throughputMap.fill(1.0);

    if(includeMechanicalVignetting  && includeOpenShutterSmearing)
        mechanicalVignettingMask.fill(1);

    double xFPmm, yFPmm;
    double angle;

//    const double refAnglePolarizationRadians = deg2rad(refAnglePolarization);       // Reference angle for the polarisation efficiency [radians]
//    const double acosPolarizationEfficiency = acos(polarizationEfficiency);

//    const double refAngleQuantumEfficiencyRadians = deg2rad(refAngleQE);     // Reference angle for the quantum efficiency [radians]
//    const double acosQuantumEfficiency = acos(relativeRefEfficiencyQE);        // Relative efficiency due to the angle dependency of the QE at the reference angle

    if (includeNaturalVignetting || includeMechanicalVignetting || includePolarization || includeQuantumEfficiency)
    {
        // Loop over all pixels in the pixel map

        for (unsigned int row = 0; row < numRowsPixelMap; row++)
        {
            for (unsigned int column = 0; column < numColumnsPixelMap; column++)
            {
                // Pixel coordinates (in the detector) -> focal-plane coordinates

                tie(xFPmm, yFPmm) = pixelToFocalPlaneCoordinates(row + subFieldZeroPointRow, column + subFieldZeroPointColumn);

                // Angular distance [radians] of the pixel from the optical axis

                angle = camera.getGnomonicRadialDistanceFromOpticalAxis(xFPmm, yFPmm);

                // Mechanical vignetting

                if(includeMechanicalVignetting)
                {
                    if (angle > radiusFOV)
                    {
                        throughputMap(row, column) = 0.0;

                        if(includeOpenShutterSmearing)
                            mechanicalVignettingMask(row, column) = 0;
                    }
                }

                // Natural vignetting.
                // With a cos^2 law, the mean natural vignetting value over all pixels is 0.945. 

                if (includeNaturalVignetting) 
                    throughputMap(row, column) *= pow(cos(angle), 2);

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












/**
 * Checks whether the gain values for left- and right-hand side of the detector are not
 * too far apart, according to the specified allowed difference.  In case the gain values
 * are too far apart, a warning message is shown to the user.
 */
void Detector::checkGain()
{
    double allowedDifference = min(refValueGainLeft, refValueGainRight) * gainAllowedDifference / 100.0;

    if(abs(refValueGainLeft - refValueGainRight) > allowedDifference)
    {
        Log.warning("Detector: Difference in gain between the left- and right-hand side of the detector too large.");
    }
}










/**
 *\brief Calculates the coefficients a^X_ij for the brighter-fatter effect,
 *       following the method proposed in Sect. 6.1 in Guyonnet et al. 2015.
 *
 *       These parameters will be the same for all pixels (0, 0) and will only
 *       be calculated for the pixels (i, j) that are within a window centred
 *       at pixel (0, 0).  This can be done because the influence of pixels (i, j)
 *       rapidly decreases with distance from pixel (0, 0).
 */
void Detector::generateGuyonnetCoefficients()
{
    Log.info("Detector: generating Guyonnet BFE coefficients.");

    // For each pixel (0, 0), we only account for the influence of pixels (i, j)
    // that are within the given range (to evaluate Eq. (11) in Guyonnet et al.).  This
    // range is defined by a window with dimensions 2 * rangeBFE + 1.

    int windowDim = 2 * rangeBFE + 1;

    // Consider the 4 directly adjacent pixels of pixel (0, 0)
    // X = {(0, 1), (0, -1), (1, 0), (-1, 0)}

    unsigned constexpr int numNeighbors = 4;
    int neighbors[numNeighbors][2] = { { 0, 1 }, { 0, -1 }, { 1, 0 }, { -1, 0 } };

    // Calculate the coefficients a^X_ij in Eq. (11) using the Eqs. in Sect. 6.1
    // in Guyonnet et al. 2015

    guyonnetCoefficients = arma::zeros<arma::Cube<float>>(windowDim, windowDim, numNeighbors);    // a^X_ij
    double rowij, columnij, r, cosTheta, f;

    // Loop over the boundaries with neighbours X

    for (unsigned int neighbor = 0; neighbor < numNeighbors; neighbor++) {

        // Loop over all "influential" pixels that are within the window

        for (int row = 0; row < windowDim; row++)
        {
            for (int column = 0; column < windowDim; column++)
            {
                rowij = (row - rangeBFE) - 0.5 * neighbors[neighbor][0];
                columnij = (column - rangeBFE) - 0.5 * neighbors[neighbor][1];

                // Distance from the source charge (Q_ij) to the considered boundary

                r = sqrt(pow(rowij, 2) + pow(columnij, 2));

                // Cosine of the angle between the source-boundary vector and the
                // normal to the boundary

                cosTheta = (rowij * neighbors[neighbor][0]
                        + columnij * neighbors[neighbor][1]) / r;

                // Eq. (18) in Guyonnet et al. 2015

                f = p0BFE * Mathematics::expint(p1BFE * r);

                guyonnetCoefficients(row, column, neighbor) = f * cosTheta;
            }
        }
    }

    // Enforcing Eq. (7) in Guyonnet et al.

    for(unsigned int neighbor = 0; neighbor < numNeighbors; neighbor++)
    {
        double sumForNeighbor = arma::accu(guyonnetCoefficients.slice(neighbor)) / pow(windowDim, 2);

        for (unsigned int row = 0; row < windowDim; row++)
        {

            for(unsigned int column = 0; column < windowDim; column++)
            {
                guyonnetCoefficients(row, column, neighbor) -= sumForNeighbor;
            }
        }

        // Accounting for the fact that (p0, p1) holds for the reference flux -> done in applyBFE()
    }
}













/**
 * \brief Verify if a point with given planar focal plane coordinates is in the subfield
 * 
 * \param xFP    Planar focal plane x-coordinate in the FP reference frame [mm]
 * \param yFP    Planar focal plane y-coordinate in the FP reference frame [mm]
 * 
 * \return true if the point is in the subfield on the CCD, false otherwise.
 */

bool Detector::isInSubfield(double xFP, double yFP)
{
    double row, column;
    tie(row, column) = focalPlaneToPixelCoordinates(xFP, yFP);

    // Check wether these pixel coordinates falls on the subfield

    return    (column >= subFieldZeroPointColumn) && (column < subFieldZeroPointColumn + numColumnsPixelMap)
           && (row    >= subFieldZeroPointRow)    && (row < subFieldZeroPointRow + numRowsPixelMap);
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

bool Detector::isInPixelMap(double row, double column)
{
    return (column >= 0) && (row >= 0) && (column < numColumnsPixelMap) && (row < numRowsPixelMap);
}











/**
 * \brief Apply throughput efficiency. This is the combined effect of:
 *          - vignetting (brightness attenuation towards the edges of the FOV);
 *          - particulate contamination;
 *          - molecular contamination;
 *          - quantum efficiency.
 */

void Detector::applyThroughputEfficiency()
{

    // Generate the throughput map which includes vignetting, polarisation, 
    // particulate & molecular contamination, and quantum efficiency. The 
    // throughput may be time dependent and therefore needs to regenerated
    // every exposure. We assume that during the jittering in Camera::exposeDetector()
    // the time dependent parameters have been updated so that their current
    // value corresponds to the very last jitter step, which are the same values
    // that we need here.
    // We also assume that the throughput does not change significantly _within_ 
    // each exposure, so that we don't need to include the throughput generation
    // in the jitter loop.

    generateThroughputMap();

    // Element-wise multiplication with the throughput map
    // Beware of Armadillo's quirky notation...
    
    pixelMap = pixelMap % throughputMap;
    
}





/**
 *\brief Adds the dark signal to the pixel map.  This follows a normal distribution,
 *       centered around the dark current and with the DSNU (percentage of the dark signal)
 *       as standard deviation.  The noise of the simulated dark signal is of the order
 *       of the sqrt of the dark signal (simulated with a normal distribution).
 *
 * \param exposureTime: Exposure time [s].
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre Bias register and smearing maps filled with zeroes.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map.
 */
void Detector::addDarkSignal(float exposureTime)
{
    double darkSignal;

    // Add dark signal to the pixel map

    // When is dark current accumulated for the pixel map?
	// 	-  exposure + readout

    const double darkCurrentOverDeltaTemp = darkCurrentStability * (getTemperature() - nominalOperatingTemperature);
    double darkSignalRef = (darkCurrent + darkCurrentOverDeltaTemp) * (exposureTime + readoutTimeBeforeNextExposure + readoutTimeDuringNextExposure);
    darkSignalDistribution = normal_distribution<double>(darkSignalRef, darkSignalRef * dsnu / 100.0);


    for(unsigned int row = 0; row < numRowsPixelMap; row++)
    {
        for(unsigned int column = 0; column < numColumnsPixelMap; column++)
        {
            // Take DSNU into account

            darkSignal = darkSignalDistribution(darkSignalGenerator);

            // Add the noise

            darkNoiseDistribution = normal_distribution<double>(darkSignal, sqrt(darkSignal));
            darkSignal = darkNoiseDistribution(darkNoiseGenerator);


            pixelMap(row, column) += darkSignal;
        }
    }

    // Add dark signal to the smearing map

    // When is dark current accumulated for the smearing map?
    // 	- normal camera, nominal mode: whole readout
    //  - fast camera, nominal mode: readout during the next exposure
    //  - partial readout: no smearing map

	darkSignalRef = darkCurrent
			* (isFastCamera ? readoutTimeDuringNextExposure : readoutTimeBeforeNextExposure);
    darkSignalDistribution = normal_distribution<double>(darkSignalRef, darkSignalRef * dsnu / 100.0);

    for(unsigned int row = 0; row < numRowsSmearingMap; row++)
    {
        for(unsigned int column = 0; column < numColumnsPixelMap; column++)
        {
            // Take DSNU into account

            darkSignal = darkSignalDistribution(darkSignalGenerator);

            // Add the noise

            darkNoiseDistribution = normal_distribution<double>(darkSignal, sqrt(darkSignal));
            darkSignal = darkNoiseDistribution(darkNoiseGenerator);

            smearingMap(row, column) += darkSignal;
        }
    }
}










/**
 * \brief: Reads out the detector and apply the following effects:
 *          - photon noise
 *          - full-well saturation (i.e. blooming)
 *          - CTE
 *          - open-shutter smearing
 *          - readout noise
 *          - quantisation:
 *              - gain
 *              - electronic offset (i.e. bias)
 *              - digital saturation
 *
 * \param exposureTime: Exposure time [s].
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre Pixel unit in the smearing map: [electrons].
 * \pre Bias register map filled with zeroes.
 *
 * \post Pixel unit in the pixel, bias register, and smearing maps: [ADU].
 */
void Detector::readOut(float exposureTime)
{

    // Add cosmic hits
    // Pixel units before: [electrons]
    // Pixel units after: [electrons]
   
    if(includeCosmicsInSubField | includeCosmicsInBiasMap | includeCosmicsInSmearingMap)
    {
        Log.debug("Detector: including cosmic hits.");
        addCosmics(exposureTime);
    }
    else
    {
        Log.debug("Detector: no cosmic hits included.");
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

    // Simulate the effects of the Charge Transfer Inefficiency (CTI). When the
    // CCD is read out, row after row, a part of the charge is always left behind
    // which then dribbles into the trailing pixels. This causes each star to have
    // a small "tail". Only visible when the CTI = 1 - CTE is poor.
    // Pixel units before: [electrons]
    // Pixel units after: [electrons]

    if (includeCTIeffects)
    {
        Log.debug("Detector: applying charge transfer inefficiency.");
        applyCTI();
    }
    else
    {
        Log.debug("Detector: no charge transfer inefficiency applied.");
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

    if(isFastCamera && frontEndElectronics->getIncludeOverAndUnderShoot())
    {
        Log.debug("Detector: adding (F-)FEE over-/undershoot");
        applyOverAndUnderShoot();
    }
    else{
        Log.debug("Detector: (F-)FEE over-/undershoot not applied: " + to_string(isFastCamera) + " " + to_string(frontEndElectronics->getIncludeOverAndUnderShoot()));
    }

    // Each time the amplifier reads out a pixel, a tiny bit of noise is added.
    // Add the readout noise.
    // Pixel units before: [electrons]
    // Pixel units after: [electrons]

    if (includeReadoutNoise)
    { 
        Log.debug("Detector: adding readout noise of CCD and FEE.");
        addReadoutNoise();
    }
    else
    {
        Log.debug("Detector: no readout noise added.");
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

    //  Apply quantisation. This consists of: 
    //         - applying FEE and CCD gain (converting from electrons to ADU)
    //         - adding the electronic offset
    //         - applying digital saturation
    // Pixel units before: [electrons]
    // Pixel units after: [ADU]
 
    
    if(includeQuantisation)
    {
        Log.debug("Detector: applying quantisation.");
        applyQuantisation();
    }
    else
    {
        Log.debug("Detector: no quantisation applied.");
    }
}









/**
 * \brief Adds the Brighter-Fatter Effect (BFE) to the pixel map, following the method
 *        proposed by Guyonnet et al. 2015 (https://arxiv.org/abs/1501.01577).
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre No bias register or smearing maps.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post No bias register or smearing maps.
 */
void Detector::applyBFE()
{
    // For each pixel (0, 0), we only account for the influence of pixels (i, j)
    // that are within the given range (to evaluate Eq. (11) in Guyonnet et al.)

    int windowDim = 2 * rangeBFE + 1;

    // Consider the 4 directly adjacent pixels of pixel (0, 0)
    // X = {(0, 1), (0, -1), (1, 0), (-1, 0)}

    unsigned constexpr int numNeighbors = 4;
    int neighbors[numNeighbors][2] = {{0, 1}, {0, -1}, {1, 0}, {-1, 0}};

    // Charges Q_0,0, Q_X, and Q_i,j in Eq. (11) in Guyonnet et al. 2015

    double charge00, chargeX, bfeX00, bfe00;//, chargeij;

    // ∂Q_0,0 in Eq. (11) in Guyonnet et al. 2015

    arma::Mat<float> bfe = arma::zeros<arma::Mat<float>>(numRowsPixelMap, numColumnsPixelMap);

    // Loop over all pixels in the pixel map that are affected by the BFE

    for(unsigned int row = rangeBFE; row < numRowsPixelMap - rangeBFE; row++)
    {
        for(unsigned int column = rangeBFE; column < numColumnsPixelMap - rangeBFE; column++)
        {
            charge00 = pixelMap(row, column);    // Q_0,0

            bfe00 = 0;

            for(unsigned int neighbor = 0; neighbor < numNeighbors; neighbor++)
            {
                chargeX = pixelMap(row + neighbors[neighbor][0], column + neighbors[neighbor][1]);    // Q_X

                // i = row - range,..., row + range
                // j = column - range, column + range

//                chargeij = pixelMap(row, column); //pixelMap(arma::span(row - rangeBFE, row + rangeBFE), arma::span(column - rangeBFE, column + rangeBFE));

                // Eq. (11) in Guyonnet et al. 2015 (within dividing by 4)
                // a^X_i,j * Q_i,j * (Q_0,0 + Q_X)
                // Accounting for the fact that (p0, p1) holds for the reference flux

                bfeX00 = arma::accu(guyonnetCoefficients.slice(neighbor) % pixelMap(arma::span(row - rangeBFE, row + rangeBFE), arma::span(column - rangeBFE, column + rangeBFE)))  / (refFluxBFE / 2.0) * (charge00 + chargeX);
                bfe00 += bfeX00;

//                for(unsigned int i = 0; i < windowDim; i++)
//                {
//                    for(unsigned int j = 0; j < windowDim; j++)
//                    {
//                        chargeij = pixelMap(row + (i - rangeBFE) , column + (j - rangeBFE));        // Q_i,j
//
//                        // Eq. (11) in Guyonnet et al. 2015 (within dividing by 4)
//                        // a^X_i,j * Q_i,j * (Q_0,0 + Q_X)
//
//                        bfe(row, column) += guyonnetCoefficients(i, j, neighbor) * chargeij / (refFluxBFE / 2.0) * (charge00 + chargeX);
//
//                        Log.info("BFE: " + to_string(guyonnetCoefficients(i, j, neighbor)) + " " + to_string(chargeij) + " " + to_string(charge00) + " " + to_string(chargeX));
//                    }
//                }
            }

            bfe(row, column) = bfe00;
        }
    }

    // Dividing by 4 in Eq. (11) in Guyonnet et al. 2015

    bfe /= 4.0;

    pixelMap += bfe;
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

    Log.debug("Adding photon noise to pixel map");

    for (unsigned int row = 0; row < numRowsPixelMap; row++)
    {
        for (unsigned int column = 0; column < numColumnsPixelMap; column++)
        {
            photonNoiseDistribution = poisson_distribution<long>(pixelMap(row, column));
            pixelMap(row, column) = photonNoiseDistribution(photonNoiseGenerator);
        }
    }

    // Add photon noise to the smearing map

    Log.debug("Adding photon noise to smearing map");

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
 * \brief: Add cosmic hits to the pixel, bias register, and smearing map.
 *         - The number of cosmic hits is determined by a random sample from a Poisson
 *           distribution with the configured mean cosmic hit rate, exposure time,
 *           size of the map (number of rows and columns [pixels]), and pixel size.
 *         - The entry points are uniformly distributed over the maps.
 *         - The entry angles are uniformly distributed over the [0, 2π] interval.
 *         - The length of the trails is uniformly distributed over the given interval.
 *         - The total number of electrons in the trail is uniformly distributed over
 *           the given interval.
 * 
 * This function is a wrapper function for Detector::addCosmics(exposureTime, map, numRows, numColumns)
 * that does the actual work but which is not called directly.
 *
 * \param exposureTime: Exposure time [s].
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre Pixel unit in the smearing map: [electrons].
 * \pre No bias register map.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post Pixel unit in the bias register map: [electrons].
 */
void Detector::addCosmics(float exposureTime)
{
	cosmicHitRateDistribution     = poisson_distribution<long>(cosmicHitRate);                                       // [hits/cm^2/s]
    cosmicEntryColumnDistribution = uniform_real_distribution<double>(0, numColumnsPixelMap - 1);                    // [pixels]
    cosmicEntryAngleDistribution  = uniform_real_distribution<double>(0, 2 * PI);                                    // [radians]
    cosmicTrailLengthDistribution = uniform_real_distribution<double>(cosmicTrailLength[0], cosmicTrailLength[1]);   // [pixels]
    cosmicIntensityDistribution   = uniform_real_distribution<double>(cosmicIntensity[0], cosmicIntensity[1]);       // [e-/hit]

    // Cosmics in the subfield

    if (includeCosmicsInSubField)
    {
        Log.debug("Detector: adding cosmic hits to the sub-field");
        addCosmics(exposureTime + readoutTimeBeforeNextExposure + readoutTimeDuringNextExposure, pixelMap, numRowsPixelMap, numColumnsPixelMap, "image area");
    }

    // Cosmics in the over-scan

    if (includeCosmicsInSmearingMap)
    {
        Log.debug("Detector: adding cosmic hits to smearing map");

        if(isFastCamera)
        {
        	addCosmics(readoutTimeDuringNextExposure, smearingMap, numRowsSmearingMap, numColumnsPixelMap, "smearing map");
        } else
        {
        	addCosmics(readoutTimeBeforeNextExposure, smearingMap, numRowsSmearingMap, numColumnsPixelMap, "smearing map");
        }
    }

    // Cosmics in the pre-scan
    // This is a special case because the rows of the prescan are all virtual. 
    // The following is only approximative. 

    if (includeCosmicsInBiasMap)
    {
        Log.debug("Detector: adding cosmic hits to bias map");
        cosmicTrailLengthDistribution = uniform_real_distribution<double>(0.0, 1.e-6);    // Only hot pixels, no trails
        const double biasMapRowLifeTime = (numColumns / 2 + numColumnsBiasMap) * serialTransferTime + parallelTransferTime;
        addCosmics(biasMapRowLifeTime, biasMapLeft, numRowsBiasMap, numColumnsBiasMap, "bias map (left half)");
        addCosmics(biasMapRowLifeTime, biasMapRight, numRowsBiasMap, numColumnsBiasMap, "bias map (right half)");
    }
}











/**
 * \brief: Add cosmic hits to the given map.
 *         - The number of cosmic hits is determined by a random sample from a Poisson
 *           distribution with the configured mean cosmic hit rate, exposure time,
 *           size of the map (number of rows and columns [pixels]), and pixel size.
 *         - The entry points are uniformly distributed over the maps.
 *         - The entry angles are uniformly distributed over the [0, 2π] interval.
 *         - The length of the trails is uniformly distributed over the given interval.
 *         - The total number of electrons in the trail is uniformly distributed over
 *           the given interval.
 *
 * This function is not called directly in Detector, but only through the method
 *      Detector::addCosmics(exposureTime)
 * 
 * \param exposureTime: amount of time exposed to cosmic particle influx [s].
 * \param map: Map affected by cosmics [e-].  Either the pixel, bias register, or
 *             smearing map.
 * \param numRows: Number of rows in the map [pixels].
 * \param numColumns: Number of columns in the map [pixels].
 * \param area: Name of the area to which the cosmics are added ("image area", "smearing map", "bias map").
 *
 * \pre Pixel unit in the pixel map: [electrons].
 * \pre Pixel unit in the smearing map: [electrons].
 * \pre No bias register map.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post Pixel unit in the bias register map: [electrons].
 */
void Detector::addCosmics(float exposureTime, arma::Mat<float> &map, int numRows, int numColumns, string area)
{
    // Characteristics of an individual trail

    double entryRow, entryColumn, entryAngle, trailLength, intensity, sigma;

    // Each trail consists of 200 "trail points"
    /// (equally spaced points over the trail length)

    int numTrailPoints = 200;
    arma::vec trailRows, trailColumns, trailWeights;    // All trail points: row, column, and weight
    int trailRow, trailColumn;                        // Individual trail point: row and column

    cosmicEntryRowDistribution = uniform_real_distribution<double>(0, numRows - 1);

    // Number of cosmic hits
    // - cosmic hit rate [events / cm^2 / s]
    // - exposure time [s]
    // - dimensions [pixels] -> [micron] -> [cm]
    //
	// To make sure the number of cosmics is as expected when considering a large number
    // of exposures, we have to account for the decimal part of the number of cosmic hits
    // per exposure too, instead of rounding down this value. E.g. if you have 1.5 cosmics
    // per exposure, you want 1500 cosmics for 1000 exposures, instead of 1000 cosmics.
    // The solution is to add as many cosmics as denoted by the integer part, and add one
    // extra cosmic, with a chance equal to the decimal part.  This is handled by the
    // uniform random distribution in [0,1]:
    // 	- random number < decimal part => add one extra cosmic
    //  - random number >= decimal part => don't add an extra cosmic

    double numCosmicHitsAsDouble = cosmicHitRateDistribution(
			cosmicHitRateGenerator) * exposureTime
			* (numRows * pixelSize / 10000.0)
			* (numColumns * pixelSize / 10000.0);
    double decimalPartNumCosmicHits = numCosmicHitsAsDouble - (int) numCosmicHitsAsDouble;

    // Round down the number of cosmic hits to an integer, possibly zero.

    int numCosmicHits = (int) numCosmicHitsAsDouble;

  	// Add 1 cosmic with a chance equal to the decimal part.

    if (decimalNumCosmicHitsDistribution(decimalNumCosmicHitsGenerator) < decimalPartNumCosmicHits)
    {
    	numCosmicHits += 1;
    }

    Log.debug("Detector: number of cosmic hits for the " + area + ": "  + to_string(numCosmicHits));
    if (numCosmicHits == 0) return;
    
    double meanEntryAngle = 0.0;
    double meanTrailLength = 0.0;
    double meanIntensity = 0.0;

    for (unsigned int cosmicHit = 0; cosmicHit < numCosmicHits; cosmicHit++)
    {
        entryRow    = cosmicEntryRowDistribution(cosmicEntryRowGenerator);          // Entry row [pixels] (uniform distribution over the rows of the sub-fields)
        entryColumn = cosmicEntryColumnDistribution(cosmicEntryColumnGenerator);    // Entry column [pixels] (uniform distribution over the columns of the sub-field)
        entryAngle  = cosmicEntryAngleDistribution(cosmicEntryAngleGenerator);      // Entry angle [radians] (uniform distribution between 0 and 2π)
        trailLength = cosmicTrailLengthDistribution(cosmicTrailLengthGenerator);    // Trail length [pixels] (uniform distribution over interval)
        intensity   = cosmicIntensityDistribution(cosmicIntensityGenerator);        // Number of e- in cosmic hit [e-] (uniform distribution over interval)

        meanEntryAngle += entryAngle;
        meanTrailLength += trailLength;
        meanIntensity += intensity;

        double trailStep = trailLength / numTrailPoints;                            // Distance between two "trail points" (for the current trail)

        trailWeights = arma::linspace(0, trailLength, numTrailPoints);
        trailRows = entryRow + trailWeights * sin(entryAngle);
        trailColumns = entryColumn + trailWeights* cos(entryAngle);

        // Apply the decay function

        sigma = arma::max(trailWeights) / 3.0;
        trailWeights = arma::exp(-arma::pow(trailWeights, 2) / (2 * pow(sigma, 2)));
        trailWeights /= arma::sum(trailWeights);

        // Add the flux coming from the cosmic hit to all its trail points in the pixel map

        for (unsigned int index = 0; index < numTrailPoints; index++)
        {
            trailRow = int(floor(trailRows(index)));
            trailColumn = int(floor(trailColumns(index)));

            if ((trailRow >= 0) && (trailRow < numRows) && (trailColumn >= 0) && (trailColumn < numColumns))
            {
                map(trailRow, trailColumn) += (trailWeights(index) * intensity);
            }
        }
    }

    meanEntryAngle /= numCosmicHits;
    meanTrailLength /= numCosmicHits;
    meanIntensity /= numCosmicHits;

    Log.info("Detector: mean cosmic entry angle (" + area + "): " + to_string(meanEntryAngle) + " rad");
    Log.info("Detector: mean cosmic trail length (" + area + "): " + to_string(meanTrailLength) + " pix");
    Log.info("Detector: mean cosmic hit intensity (" + area + "): " + to_string(meanIntensity) + " e-");
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
 * \pre No bias register map, unless cosmics have been added to it.
 * \pre Full-well saturation limit expressed in [electrons].
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Effect of full-well saturation (i.e. blooming) applied to the pixel map.
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map, unless cosmics have been added to it.
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
 * \pre No bias register map, unless cosmics have been added to it.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map, unless cosmics have been added to it.
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
 * \pre No bias register map, unless cosmics have been added to it.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map, unless cosmics have been added to it.
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

    arma::Row<float> readout;   // Readout strip

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
                const double cteFactor = pow(meanCte, index + 1) * pow(cti, row + subFieldZeroPointRow - index);

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
 * \pre No bias register map, unless cosmics have been added to it.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map, unless cosmics have been added to it.
 */
 
void Detector::applyShort2013CTImodel()
{
    // Compute the maximum geometrical volume that electrons can occupy within a pixel.
    // I.e. the volume of the electron cloud when the capacity of the full well is maxed out.
    // E.g if the pixel size is 18 micron, then the volume is 18e-6 * 18e-6 * 1.e-6 / 2.0

    const double maxVolumePerPixel = pixelSize * pixelSize * 1.e-18  / 2.0;                                  // Vg [m^3]

    // Compute the time it takes to transfer 1 row during readout

    const double chargeTransferTime = parallelTransferTime;		// t[s]
//    const double chargeTransferTime = readoutTime / numRows;                                                 // t [s]

    // Compute the thermal velocity of the electrons in the silicon

    const double effectiveElectronMass = 0.5 * Constants::FREEELECTRONMASS;                                  // me [kg]
    const double thermalVelocity = sqrt(3.0 * Constants::KBOLTZMANN * temperature / effectiveElectronMass);  // vt [m/s]

    // Arrays to keep track of the number of occupied traps in a column

    arma::Mat<float> numberOfOccupiedTraps = arma::zeros<arma::Mat<float>>(numTrapSpecies, numColumnsPixelMap);	// No

    // Arrays to keep track of the captured and released electrons, for each column in a particular row.

    arma::Row<float> numberOfCapturedElectrons(numColumnsPixelMap);		// Nc
    arma::Row<float> numberOfReleasedElectrons(numColumnsPixelMap);		// Nr

    arma::Row<float> alpha(numTrapSpecies, arma::fill::zeros);

    // Eq. (23) of Short et al. 2013

    for (int k = 0; k < numTrapSpecies; k++)
    {
    		alpha(k) = chargeTransferTime * trapCaptureCrossSection[k] * thermalVelocity * pow(fullWellSaturationLimit, beta) / (2.0 * maxVolumePerPixel);
    }

    // Loop over all rows of the pixelMap, and over all trap species.
    // For each row, the computations are done for all columns simultaneously.

    for (int rowNumber = 0; rowNumber < numRowsPixelMap; rowNumber++)
    {

        for (int k = 0; k < numTrapSpecies; k++)
        {
            // Compute the number of electrons captured in a trap, according to Eq. (22)-(23) of Short et al. (2013).
            // Note that Armadillo uses % for elementwise multiplication.

            const double gamma = 2 * trapDensity[k] * maxVolumePerPixel / pow(fullWellSaturationLimit, beta) * (subFieldZeroPointRow + rowNumber + 1);	// +1 as row = 0 also has to be transferred once

            numberOfCapturedElectrons =   (gamma * arma::pow(pixelMap.row(rowNumber), beta) - numberOfOccupiedTraps.row(k)) \
                                        / (gamma * arma::pow(pixelMap.row(rowNumber), beta-1) + 1)                          \
                                        % (1 - arma::exp(-alpha(k) * arma::pow(pixelMap.row(rowNumber), 1-beta)));

            // Captured electron numbers can't be negative, so clip negative value to zero.

            arma::Col<arma::uword> isNegative = arma::find(numberOfCapturedElectrons < 0.0);
            numberOfCapturedElectrons(isNegative).zeros();
            numberOfCapturedElectrons.elem(find_nonfinite(numberOfCapturedElectrons)).zeros();

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
 * \pre No bias register map, unless cosmics have been added to it.
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post No bias register map, unless cosmics have been added to it.
 */
void Detector::applyOpenShutterSmearing(float exposureTime)
{
    if (includeMechanicalVignetting)
    {
        // The mask indicating which pixels in the sub-field are within the FOV and which
        // ones not, has already been created upon construction of the throughput map.

        double rowFOV, xFPmmFirstExposed, yFPmmFirstExposed, xFPmmTop, yFPmmTop, angleFirstExposed, angleTop;

        // Loop over all columns in the sub-field

        for(unsigned int column = 0; column < numColumnsPixelMap; column++) //[pixels in sub-field]
        {
            // Intersection of the current column of the detector with the circle representing the FOV:
            // (rowFOV, column).  If no intersection can be found, this is NaN.
            rowFOV = getRowEdgeFOV(column);

            if(isnan(rowFOV))
            {
                // All exposed rows are within the FOV

                numExposedRowsInFOV(column) = numRows - firstRowExposed;
            }
            else
            {
                tie(xFPmmFirstExposed, yFPmmFirstExposed) = pixelToFocalPlaneCoordinates(firstRowExposed, column + subFieldZeroPointColumn);    // detector [pixels] -> focal plane [mm]
                angleFirstExposed = camera.getGnomonicRadialDistanceFromOpticalAxis(xFPmmFirstExposed, yFPmmFirstExposed);

                tie(xFPmmTop, yFPmmTop) = pixelToFocalPlaneCoordinates(numRows - 1, column + subFieldZeroPointColumn);
                angleTop = camera.getGnomonicRadialDistanceFromOpticalAxis(xFPmmTop, yFPmmTop);

                if (angleFirstExposed > radiusFOV)
                {
                    if (angleTop > radiusFOV)
                        numExposedRowsInFOV(column) = 0;
                    
                    else
                        numExposedRowsInFOV(column) = numRows - (int) rowFOV;
                }

                else
                {
                    if (angleTop > radiusFOV)
                        numExposedRowsInFOV(column) = (int) rowFOV - firstRowExposed + 1;
                    else
                        numExposedRowsInFOV(column) = numRows - firstRowExposed;
                }
            }
        }
    }



    // Average out the fluxes in the pixel map per column (of the whole CCD) and make sure it is
    // scaled with the readout time (during which the detector is susceptible to open-shutter smearing) 
    // instead of with the exposure time:
    // - rows in the sub-field (in the exposed part of the detector and in the FOV): use actual fluxes (accumulated during exposure)
    // - rows outside the sub-field (in the exposed part of the detector and in the FOV): use total sky background
    // - rows outside the FOV and/or not exposed are not considered (as these rows are shielded off against incoming radiation)

    arma::Row<float> openShutterSmearing = arma::sum(pixelMap % mechanicalVignettingMask, 0);   // Flux in the exposed part of the sub-field that is inside the FOV (per column)
    arma::Row<int> numExposedSubFieldRowsInFOV = arma::sum(mechanicalVignettingMask, 0);        // Number of pixels in the exposed part of the sub-field that are inside the FOV (per column)

    arma::Row<int> numExposedNonSubFieldRowShifts = numExposedRowsInFOV - numExposedSubFieldRowsInFOV + numRowsSmearingMap; // Number of pixels in the exposed part of the detector that do not reside in the sub-field + parallel over-scan
    arma::Row<float> openShutterSmearingOutsideSubField = arma::conv_to<arma::Row<float>>::from(numExposedNonSubFieldRowShifts) * camera.getTotalSkyBackground();

   // Apply all throughput efficiencies (these have not been applied to the total sky background yet)
   // Note that mechanical vignetting has already been taken into account (if applicable)

    if(includeNaturalVignetting)
        openShutterSmearingOutsideSubField *= expectedValueNaturalVignetting;

    if(includePolarization)
        openShutterSmearingOutsideSubField *= expectedValuePolarization;

    if(includeParticulateContamination)
        openShutterSmearingOutsideSubField *= particulateContaminationEfficiency;

    if(includeMolecularContamination)
        openShutterSmearingOutsideSubField *= molecularContaminationEfficiency;

    // Sect. 4.2.4.5 of PLATO-DLR-PL-RP-001:
    // The expected value E_ang is then the mean over all pixels and results in a value of 0.993

    if(includeQuantumEfficiency)
        openShutterSmearingOutsideSubField *= meanQE;

    openShutterSmearing += openShutterSmearingOutsideSubField;

    // Fast camera: lower half of the CCD is shielded off -> no contribution to open-shutter smearing here

    arma::Row<float> openShutterSmearingTime = arma::conv_to<arma::Row<float>>::from(numExposedRowsInFOV + numRowsSmearingMap) / 
        (numRows - firstRowExposed + numRowsSmearingMap) * readoutTimeBeforeNextExposure;

    arma::Row<float> factor = openShutterSmearingTime % arma::pow(arma::conv_to<arma::Row<float>>::from(numExposedRowsInFOV), -1) / exposureTime;

    openShutterSmearing %= factor;

    // Add the effect of the open-shutter smearing to the pixel map

    for (unsigned int row = 0; row < numRowsPixelMap; row++)
        pixelMap(row, arma::span::all) += openShutterSmearing;

    // Add the effect of the open-shutter smearing to the smearing map

    for (unsigned int row = 0; row < numRowsSmearingMap; row++)
        smearingMap(row, arma::span::all) += openShutterSmearing;
}










/**
 * \brief Find the intersection of the circle representing the FOV and the given column of the
 *        detector.  In case this circle does intersect with the detector, the row coordinate
 *        is returned; if not, NaN is returned.
 * 
 * \param column: Column (coordinate) [pixels], for which to find the intersection with the
 *                circle representing the FOV.
 * 
 * \return In case this circle does intersect with the detector, the row coordinate
 *         is returned; if not, NaN is returned.
 **/ 
double Detector::getRowEdgeFOV(int column)
{
    double pixelSizeMm = pixelSize / 1000.0;    // Pixel size [µm] -> [mm]

    // Quadratic equation: a * x**2 + b * x + c  = 0
    // Find intersection between circle representing the FOV and the given column

    double a = pow(pixelSizeMm, 2);
    double b = - 2 * pixelSizeMm * subFieldZeroPointRow;
    double c = pow(subFieldZeroPointRow, 2) + pow(column * pixelSizeMm - subFieldZeroPointColumn , 2) -  pow(camera.getFocalLength() * tan(radiusFOV), 2);

    // Discriminant (should be positive)

    double discriminant = pow(b, 2) - 4 * a * c;

    if (discriminant < 0)
        return nan("");
    
    double discriminantSqrt = sqrt(discriminant);

    // First solution of the equation
    // Check whether it intersects with the detector (if not, try the other solution)

    double solutionRow = (-b + discriminantSqrt) / (2 * a);

    if ((solutionRow >= 0) && (solutionRow < numRows))
        return solutionRow;

    // Second solution of the equation
    // Check whther it intersects with the detector (if not, the whole row is within the FOV)

    solutionRow = (-b - discriminantSqrt) / (2 * a);
    
    if ((solutionRow >= 0) && (solutionRow < numRows))
        return solutionRow;

    // No proper solution found (i.e. circle of FOV does not intersect
    // with the detector in the given column)

    return nan("");
}













/**
 * \brief Apply the readout noise to the pixel map, bias map, and smearing map.  The readout
 *        noise is contributed to by the detector and by the FEE.
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
 * \pre No bias register map, unless cosmics have been added to it.
 * \pre Readout noise expressed in [electrons].
 *
 * \post Pixel unit in the pixel map: [electrons].
 * \post Pixel unit in the smearing map: [electrons].
 * \post Pixel unit in the bias register map: [electrons].
 * \post Added readout noise to the bias and smearing maps.
 */

void Detector::addReadoutNoise()
{
    // Adding the readout noise in quadrature (CCD & FEE)

    const double totalReadoutNoise = sqrt(pow(readoutNoise, 2) + pow(frontEndElectronics->getReadoutNoise(), 2));

    readoutNoiseDistribution = normal_distribution<double>(0.0, totalReadoutNoise);

    // Add readout noise to the pixel map

    for (unsigned int row = 0; row < numRowsPixelMap; row++)
    {
        for (unsigned int column = 0; column < numColumnsPixelMap; column++)
        {
            pixelMap(row, column) += readoutNoiseDistribution(readoutNoiseGenerator);
        }
    }

    // Add readout noise to the bias prescan map

    for(unsigned int row = 0; row < numRowsBiasMap; row++)
    {
    	for(unsigned int column = 0; column < numColumnsBiasMap; column++)
    	{
    		biasMapLeft(row, column) += readoutNoiseDistribution(readoutNoiseGenerator);
    		biasMapRight(row, column) += readoutNoiseDistribution(readoutNoiseGenerator);
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
 * \brief: Apply quantisation to the bias register, smearing, and smearing map.  This consists of
 *         the following steps:
 *         - applying FEE and CCD gain (converting from electrons to ADU)
 *         - adding the electronic offset
 *         - forcing the ADUs to be integers
 *         - applying digital saturation
 *
 * \pre Pixel unit in the pixel, smearing, and bias register maps: [electrons].
 *
 * \post Pixel unit in the pixel, smearing, and bias register maps: [ADU].
 */
void Detector::applyQuantisation()
{
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
    biasMapLeft = arma::floor(biasMapLeft);
    biasMapRight = arma::floor(biasMapRight);
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
 * \brief: Divide the bias register, smearing, and pixel map by the detector gain.
 *         This converts these three maps from electrons to ADU.
 *
 * \pre Pixel unit in the pixel, smearing, and bias register maps: [electrons].
 *
 * \post Pixel unit in the pixel, smearing, and bias register maps: [ADU].
 */
void Detector::applyGain()
{
    Log.debug("Detector: applying gain to pixelMap, biasMap and smearingMap");

    const int lastIndexCcdLeft = numColumns / 2 - 1;
    const int lastIndexSubFieldLeft = lastIndexCcdLeft - subFieldZeroPointColumn;

    // Detector gain (left & right) [µV / e-]

    const double ccdGainOverDeltaTemp = gainStability * (getTemperature() - nominalOperatingTemperature);

    const double ccdGainLeft = refValueGainLeft + ccdGainOverDeltaTemp;
    const double ccdGainRight = refValueGainRight + ccdGainOverDeltaTemp;

    // FEE gain (left & right) [ADU / µV]

    // Combined gain (FEE & CCD) [ADU / e-]

    const double combinedGainLeft = frontEndElectronics->getGainLeftAdc(internalTime) * ccdGainLeft;
    const double combinedGainRight = frontEndElectronics->getGainRightAdc(internalTime) * ccdGainRight;

    if(lastIndexSubFieldLeft >= numColumnsPixelMap - 1)      // Left ADC only
    {
        pixelMap *= combinedGainLeft;
        smearingMap *= combinedGainLeft;
    }
    else if(lastIndexSubFieldLeft < 0)                     // Right ADC only
    {
        pixelMap *= combinedGainRight;
        smearingMap *= combinedGainRight;
    }
    else
    {
        // 0 -> lastIndexSubFieldLeft: left ADC

        pixelMap.submat(arma::span::all, arma::span(0, lastIndexSubFieldLeft)) *= combinedGainLeft;
        smearingMap.submat(arma::span::all, arma::span(0, lastIndexSubFieldLeft)) *= combinedGainLeft;

        // lastIndexSubFieldLeft + 1 -> numColumnsSubPixelMap -1: right ADC

        pixelMap.submat(arma::span::all, arma::span(lastIndexSubFieldLeft, numColumnsPixelMap - 1)) *= combinedGainRight;
        smearingMap.submat(arma::span::all, arma::span(lastIndexSubFieldLeft, numColumnsPixelMap - 1)) *= combinedGainRight;
    }

    biasMapLeft *= combinedGainLeft;
    biasMapRight *= combinedGainRight;

    Log.info("Detector: gain of left part of CCD: " + to_string(combinedGainLeft));
    Log.info("Detector: gain of right part of CCD: " + to_string(combinedGainRight));
}











/**
 * \brief: Add the electronic offset (i.e. bias level) of the FEE to the pixel map,
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
    const double offset = frontEndElectronics->getElectronicOffset(internalTime);

    Log.debug("Detector: adding a bias to pixelMap, biasMap and smearingMap (electronicOffset=" + to_string(offset)+ ")");

    // Add the electronic offset to the pixel, bias register, and smearing maps

    pixelMap += offset;
    biasMapLeft += offset;
    biasMapRight += offset;
    smearingMap += offset;
}












/**
 * \brief Apply the effect of digital saturation to the pixel map,
 *        smearing map, and bias register map. This means that the pixel values in
 *        these maps (expressed in [ADU / pixel]) are Lastped off to the digital saturation
 *        limit of the detector (also expressed in [ADU / pixel]).
 *
 * \pre Pixel unit in the pixel, smearing, and bias maps: [ADU].
 * \pre Digital saturation limit expressed in [ADU / pixel].
 *
 * \post Pixel unit in the pixel, smearing, and bias register maps: [ADU].
 */
void Detector::applyDigitalSaturation()
{
    // Last off the values in the pixel map

    pixelMap(arma::find(pixelMap > digitalSaturationLimit)).fill(digitalSaturationLimit);

    // Last off the values in the bias register map

    biasMapLeft(arma::find(biasMapLeft > digitalSaturationLimit)).fill(digitalSaturationLimit);
    biasMapRight(arma::find(biasMapRight > digitalSaturationLimit)).fill(digitalSaturationLimit);

    // Last off the values in the smearing map

    smearingMap(arma::find(smearingMap > digitalSaturationLimit)).fill(digitalSaturationLimit);
}












/**
 * \brief Applies the F-FEE over-/undershoot to the pixel map.
 * 
 * \pre Pixel unit in the pixel, smearing, and bias maps: [ADU].
 *
 * \post Pixel unit in the pixel, smearing, and bias register maps: [ADU].
 */ 
void Detector::applyOverAndUnderShoot()
{
    const unsigned int halfDectectorWidth = numColumns / 2;

    // NOTES
    // - We apply this effect row per row, because otherwise, we will have to keep too much
    //   data in memory
    // - Both detector halves must be treated independently
    // - The pixels between the bias map and the pixel map are assumed to contain the sky background

    arma::frowvec readoutRegister;   // Placeholder
    arma::frowvec difference;        // Placeholder
    arma::frowvec totalContribution; // Has to be filled with zeroes for every row and every detector half
    int lengthReadoutRegister;

    double skyBackground = camera.getTotalSkyBackground();
    if (includeNaturalVignetting)
        skyBackground *= expectedValueNaturalVignetting;
    if (includePolarization)
        skyBackground *= expectedValuePolarization;
    if (includeParticulateContamination)
        skyBackground *= particulateContaminationEfficiency;
    if (includeMolecularContamination)
        skyBackground *= molecularContaminationEfficiency;
    skyBackground *= meanQE * meanAngleDependencyQE;

    arma::frowvec skyBackgroundExtraPixels(frontEndElectronics->getOverAndUnderShootRange());
    skyBackgroundExtraPixels.fill(skyBackground);

    // At least partially on the left detector half

    if (subFieldZeroPointColumn < halfDectectorWidth)
    {
        const int firstIndexLeftHalf = subFieldZeroPointColumn;
        const int lastIndexLeftHalf = min(halfDectectorWidth - 1, subFieldZeroPointColumn + numColumnsPixelMap - 1);
        const int numCcdPixelsLeftHalf = lastIndexLeftHalf - firstIndexLeftHalf + 1;
     
        lengthReadoutRegister = numCcdPixelsLeftHalf + frontEndElectronics->getOverAndUnderShootRange();    // Pixels in sub-field on left CCD half + some extra pixels closest to the readout electronics
        readoutRegister.zeros(lengthReadoutRegister);
        readoutRegister.head(frontEndElectronics->getOverAndUnderShootRange()) = skyBackgroundExtraPixels;

        for (unsigned int row = 0; row < numRowsPixelMap; row++)
        {
            readoutRegister(arma::span(frontEndElectronics->getOverAndUnderShootRange(), lengthReadoutRegister - 1)) = pixelMap.row(row).head(numCcdPixelsLeftHalf);
            totalContribution.zeros(lengthReadoutRegister);

            difference = readoutRegister.tail(lengthReadoutRegister - 1) - readoutRegister.head(lengthReadoutRegister - 1);

            for (int deltaX = 0; deltaX < frontEndElectronics->getOverAndUnderShootRange(); deltaX++)
            {
                // Ditch last deltaX and tail(numCcdPixelsLeftHalf)

                totalContribution.tail(lengthReadoutRegister - frontEndElectronics->getOverAndUnderShootRange()) += frontEndElectronics->getOverAndUnderShootStrength() * difference.head(lengthReadoutRegister - 1 - deltaX).tail(numCcdPixelsLeftHalf) * exp(-frontEndElectronics->getOverAndUnderShootDecayRate() * pow(deltaX, frontEndElectronics->getOverAndUnderShootDecaySpeed()));
                
                // totalContribution.tail(lengthReadoutRegister - deltaX) += frontEndElectronics->getOverAndUnderShootStrength() * difference * exp(-frontEndElectronics->getOverAndUnderShootDecayRate() * pow(deltaX, frontEndElectronics->getOverAndUnderShootDecaySpeed()));
            }

            pixelMap.row(row).head(numCcdPixelsLeftHalf) += totalContribution.tail(lengthReadoutRegister - frontEndElectronics->getOverAndUnderShootRange());
        }
    }

    // At least partially on the right detector half

    if (subFieldZeroPointColumn + numColumnsPixelMap - 1 >= halfDectectorWidth)
    {
        const int firstIndexRightHalf = max(halfDectectorWidth, subFieldZeroPointColumn);
        const int lastIndexRightHalf = subFieldZeroPointColumn + numColumnsPixelMap - 1;
        const int numCcdPixelsRightHalf = lastIndexRightHalf - firstIndexRightHalf + 1;

        lengthReadoutRegister = numCcdPixelsRightHalf + frontEndElectronics->getOverAndUnderShootRange();      // Pixels in sub-field on right CCD half + some extra pixels closest to the readout electronics
        readoutRegister.zeros(lengthReadoutRegister);
        readoutRegister.tail(frontEndElectronics->getOverAndUnderShootRange()) = skyBackgroundExtraPixels;

        for (unsigned int row = 0; row < numRowsPixelMap; row++)
        {
            readoutRegister(arma::span(0, lengthReadoutRegister - frontEndElectronics->getOverAndUnderShootRange() - 1)) = pixelMap.row(row).tail(numCcdPixelsRightHalf);
            totalContribution.zeros(lengthReadoutRegister);

            difference = readoutRegister.head(lengthReadoutRegister - 1) - readoutRegister.tail(lengthReadoutRegister - 1);

            for (int deltaX = 0; deltaX < frontEndElectronics->getOverAndUnderShootRange(); deltaX++)
            {
                // Ditch first deltaX and head(numCcdPixelsRightHalf)

                totalContribution.head(lengthReadoutRegister - frontEndElectronics->getOverAndUnderShootRange()) += frontEndElectronics->getOverAndUnderShootStrength() * difference.tail(lengthReadoutRegister - 1 - deltaX).head(numCcdPixelsRightHalf) * exp(-frontEndElectronics->getOverAndUnderShootDecayRate() * pow(deltaX, frontEndElectronics->getOverAndUnderShootDecaySpeed()));
            }

            pixelMap.row(row).tail(numCcdPixelsRightHalf) += totalContribution.head(lengthReadoutRegister - frontEndElectronics->getOverAndUnderShootRange());
        }
    }
}














/**
 * \brief Compute the (x,y) coordinates [mm] in the FP reference system 
 *        given the (real-valued) pixel row and column numbers on the CCD.
 *        
 * \note  The rows correspond to the y-direction, and the columns to the x-direction.
 *        Pixel (row, col) = (0,0) starts at (yFP, xFP) = (0, 0).
 *               
 * \param row     CCD row coordinate, real-valued (e.g. 3.5)    [pix]
 * \param column  CCD column coordinate, real-valued (e.g. 8.3) [pix]
 * 
 * \return (xFP, yFP)  A pair of (x,y) coordinates in the FP reference system [mm]
 */

pair<double, double> Detector::pixelToFocalPlaneCoordinates(double row, double column)
{

    // Convert the pixel coordinates into [mm] coordinates

    const double xCCDmm = column * pixelSize / 1000.0;
    const double yCCDmm = row * pixelSize / 1000.0;

    // Convert the CCD coordinates into FP coordinates [mm]

    const double xFP = (xCCDmm - originOffsetX) * cos(orientationAngle) - (yCCDmm - originOffsetY) * sin(orientationAngle);
    const double yFP = (xCCDmm - originOffsetX) * sin(orientationAngle) + (yCCDmm - originOffsetY) * cos(orientationAngle);

    // That's it

    return make_pair(xFP, yFP);
}










/**
 * \brief Compute the (real-valued) pixel coordinates of the star on the CCD, given the 
 *        (x,y) coordinates [mm] in the FP reference system
 *
 * \note  - The rows correspond to the y-direction, and the columns to the x-direction.
 *        - Pixel (row, col) = (0,0) starts at (yFP, xFP) = (0, 0).
 *        
 * \param xFP  x-coordinate of the point in the FP reference system  [mm]
 * \param yFP  y-coordinate of the point in the FP reference system  [mm]
 * 
 * \return (row, column)  row and column pixel coordinates of the point (real-valued) [pix]
 */

pair<double, double> Detector::focalPlaneToPixelCoordinates(double xFP, double yFP)
{
    // Convert the FP coordinates into CCD coordinates [mm]

    const double xCCDmm = originOffsetX + xFP * cos(orientationAngle) + yFP * sin(orientationAngle);
    const double yCCDmm = originOffsetY - xFP * sin(orientationAngle) + yFP * cos(orientationAngle);

    // Convert the [mm] coordinates into pixel coordinates

    const double column = xCCDmm / pixelSize * 1000.0;
    const double row = yCCDmm / pixelSize * 1000.0;

    // That's it

    return make_pair(row, column);
}












/**
 * \brief  Return the focal plane coordinates of the center pixel of the subfield
 * 
 * \return (xFP, yFP)   focal plane coordinates in the FP'reference system [mm]
 */

pair<double, double> Detector::getFocalPlaneCoordinatesOfSubfieldCenter()
{
    double centerRow = subFieldZeroPointRow + numRowsPixelMap / 2.0;
    double centerCol = subFieldZeroPointColumn + numColumnsPixelMap / 2.0;

    // The columns correspond to the x-coordinate, the rows to the y-coordinate

    double xFP, yFP;
    tie(xFP, yFP) = pixelToFocalPlaneCoordinates(centerRow, centerCol);

    return make_pair(xFP, yFP);
}










/**
 * \brief Return the (X,Y) coordinates in the FP' reference frame in [mm] of the 4 corners
 *        of the subfield
 *        
 * \return (X00, Y00, X01, Y01, X11, Y11, X10, Y10)  [mm]
 *         where: (X00, Y00) are the FP coordinates of the lower left corner of the subfield
 *                (X01, Y01) are the FP coordinates of the lower right corner of the subfield
 *                (X11, Y11) are the FP coordinates of the upper right corner of the subfield
 *                (X10, Y10) are the FP coordinates of the upper left corner of the subfield
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
    hdf5File.createGroup("/BiasMapsLeft");
    hdf5File.createGroup("/BiasMapsRight");
    hdf5File.createGroup("/SmearingMaps");
    hdf5File.createGroup("/Flatfield");
    hdf5File.createGroup("/ThroughputMaps");
}












/**
 * \brief: Writes the pixel map for the HDF5 file.
 * 
 * \param exposureNr:   Sequential number of the exposure
 */
void Detector::writePixelMapsToHDF5(int exposureNr)
{
	// Compose the image name

    stringstream myStream;
    myStream << "image" << setfill('0') << setw(6) << exposureNr;
    string imageName = myStream.str();

    // Add the image to the "Images" group

    if (!includeQuantisation)
    {
        // Write the float array to HDF5

        hdf5File.writeArray("/Images", imageName, pixelMap);
    }
    else
    {
        // Write the pixel maps as 2-byte (16 bit) unsigned short integers.
        // As a safety check, first check that the extrema of the map are indeed
        // within the boundaries of such a data type.
       
        if((pixelMap.min() < 0) || (pixelMap.max() >= (1 << 16)))
        {
            throw ConfigurationException("Detector: quantisation was applied but pixel map values are not in [0, 2^16[");
        }

        // Convert the float matrix to an unsigned uint16_t matrix

        arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(pixelMap);
        hdf5File.writeArray("/Images", imageName, uintMap);
    }

    // send the imagette to the client connected via tcp connection

    if (sendImagettesToClient)
    {
        Log.info("Detector: attempt to send the imagette to the client");

        // convert the imagette to const char*

        std::string imagetteString = convertMatrixToString(&pixelMap, exposureNr);

        const char* imagetteChar = imagetteString.c_str();

        // declare the the outgoing zmq::message

        zmq::message_t imagetteMessage (strlen(imagetteChar));

        memcpy(imagetteMessage.data(), imagetteChar, strlen(imagetteChar));

        // send the imagette to the client

        imagetteSocket->send(imagetteMessage);
    }


    if (numRowsSmearingMap != 0)
    {
    	// Clear the string stream and compose the smearing map name

    	myStream.str(string());      // insert empty string
    	myStream.clear();            // clear eof bit

    	myStream << "smearingMap" << setfill('0') << setw(6) << exposureNr;
    	string smearingMapName = myStream.str();


    	// Add the smearing map to the "SmearingMaps" group

        if (!includeQuantisation)
        {
            // Write the float array to HDF5

    	    hdf5File.writeArray("/SmearingMaps", smearingMapName, smearingMap);
        }
        else
        {
            if ((smearingMap.min() < 0) || (smearingMap.max() >= (1 << 16)))
            {
                throw ConfigurationException("Detector: quantisation was applied but smearing map values are not in [0, 2^16[");
            }

            // Convert the float matrix to an unsigned uint16_t matrix

            arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(smearingMap);
            hdf5File.writeArray("/SmearingMaps", smearingMapName, uintMap);
        }
        
    }

   // Clear the string stream and compose the bias map name

    myStream.str(string());      // insert empty string
    myStream.clear();            // clear eof bit

    myStream << "biasMap" << setfill('0') << setw(6) << exposureNr;
    string biasMapName = myStream.str();

    // Add the bias map to the "BiasMaps" group

    if (!includeQuantisation)
    {
        // Write the float array to HDF5

        hdf5File.writeArray("/BiasMapsLeft", biasMapName, biasMapLeft);
        hdf5File.writeArray("/BiasMapsRight", biasMapName, biasMapRight);
    }
    else
    {
        if ((biasMapLeft.min() < 0) || (biasMapLeft.max() >= (1 << 16)))
        {
            throw ConfigurationException("Detector: quantisation was applied but pixel values in the left bias map are not in [0, 2^16[");
        }

        if ((biasMapRight.min() < 0) || (biasMapRight.max() >= (1 << 16)))
        {
            throw ConfigurationException("Detector: quantisation was applied but pixel values in the right bias map are not in [0,2^16[");
        }

        // Convert the float matrix to an unsigned uint16_t matrix

        arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(biasMapLeft);
        hdf5File.writeArray("/BiasMapsLeft", biasMapName, uintMap);

        uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(biasMapRight);
        hdf5File.writeArray("/BiasMapsRight", biasMapName, uintMap);
    }
    

    // Clear the string stream and compose the throughput map name

    myStream.str(string());      // insert empty string
    myStream.clear();            // clear eof bit

    myStream << "throughputMap" << setfill('0') << setw(6) << exposureNr;
    string throughputMapName = myStream.str();

    // Add the throughput map to the "ThroughputMaps" group

    hdf5File.writeArray("/ThroughputMaps", throughputMapName, throughputMap);
}








/**
 * \brief Returns the current temperature of the detector.
 */
double Detector::getTemperature()
{
    return temperatureGenerator.getNextTemperature(internalTime);
}












/**
 * \brief Returns the duration of the readout before the next exposure can start [s].
 */
double Detector::getReadoutTimeBeforeNextExposure()
{
	return readoutTimeBeforeNextExposure;
}




/**
 * \brief sets the imagette socket to clarify where imagettes are to be send, when finished.
 */
void Detector::setImagetteSocket(zmq::socket_t* socket)
{
    imagetteSocket = socket;
    sendImagettesToClient = true;
}


/**
 * \brief sets the window position socket to clarify where the input parameters are coming from
 */
void Detector::setWinPositionSocket(zmq::socket_t* socket)
{
    winPositionSocket = socket;
    getWinPositionFromServer = true;
}



/**
 * \brief converts the imagetteID, rows, cols and the pixelmap values to a char (seperated by a blank space) to be send to the client
 *  
 */
std::string Detector::convertMatrixToString(arma::Mat<float>* pixelMapPointer, uint imagetteCounter)
{
    int rows = pixelMapPointer->n_rows;

    int cols = pixelMapPointer->n_cols;

    // write the values seperated by a white space to the string
    std::string imagetteString = to_string(imagetteCounter) + " " + to_string(rows) + " " + to_string(cols) + " ";

    // write every value of the pixelMap to the string
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            imagetteString += to_string(int(pixelMapPointer->at(i, j))) + " ";
        }
    }

    return imagetteString;
}


/**
 * \brief: checks whether there is a message to change the window postion
 *         the first time the application reaches this point it waits for a message from the server
 *         after that, it only checks, whether there is news and carries on with the old data if not
 *  
 */
void Detector::setWinPosition()
{

    // send a handshake message to the position server and wait for a response before the simulation starts

    if(firstExposure)
    {
        std::string messageString = "";

        zmq::message_t message(messageString.length());

        const char *cMessage = messageString.c_str(); 

        memcpy (message.data (), cMessage, messageString.length());


        winPositionSocket->send(message);


        Log.info("DetectorWithMappedPSF: Wait for Message from window position server");

        // wait for message from winListServer

        zmq::message_t reply;

        std::string replyString;

        winPositionSocket->recv(&reply);

        replyString = std::string(static_cast<char*>(reply.data()), reply.size());

        // do a sanity check of the received message and set the window list

        if(checkWinPositionMessage(replyString))
        {
            Log.info("DetectorWithMappedPSF: got first win position message");

            // if the window position is set the simulation can start

            firstExposure = false;

            // set the socket to no longer wait for input from server

            uint timeOut = 0;

            winPositionSocket->setsockopt(ZMQ_RCVTIMEO, &timeOut, sizeof(timeOut));

        }

    }

    Log.info("DetectorWithMappedPSF: Check whether there is a new Message from window position server");

    // after the simulation is started by the reply of the server (and the corresponding setting of the window position)
    // the platoSim instance should check whether there are new commands for each new imagette creation

    // define a message
    zmq::message_t winPositionMessage;
  
    // check whether there is a new message
    if (winPositionSocket->recv(&winPositionMessage))
    {
        Log.info("DetectorWithMappedPSF: Got Message from window position server");

        std::string messageString = std::string(static_cast<char*>(winPositionMessage.data()), winPositionMessage.size());

        bool newPositionSet = checkWinPositionMessage(messageString);                
    }
    else
    {
        Log.info("Detector: no message received");
    }

}

bool Detector::checkWinPositionMessage(std::string message)
{
    
}