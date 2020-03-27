#include "DetectorWithSymmetricalMappedPSF.h"

/**
 * \brief Constructor.
 * 
 * \details Initialises the groups in the HDF5 file where the different maps (i.e. pixel map,
 *          bias register map, smearing map, etc.) will be saved. 
 * 
 * The following maps are initialized to zero (partly through the base class Detector):
 *      - pixelMap 
 *      - subPixelMap
 *      - biasMap
 *      - smearingMap
 *      - flatfieldMap
 *      - throughputMap
 *      - cteMap
 * 
 * The flatfieldMap is filled at sub-pixel level, the throughputMap and cteMap are filled at pixel level.
 * 
 * \param configParam: Configuration parameters for the detector.
 * 
 * \param hdf5file:HFD5 file to write the detector images to.
 * 
 * \param camera: Camera to which to attach the detector.
 * 
 * \param readoutTimeBeforeNextExposure Duration of the readout that takes place before the next exposure can start [s].
 */
DetectorWithSymmetricalMappedPSF::DetectorWithSymmetricalMappedPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
    : DetectorWithMappedPSF(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure)
{
    // Parse the parameters from the configuration file.

    configure(configParam);

    // Create the groups in the HDF5 file where the different maps (i.e. pixel map,
    // bias register map, smearing map, etc.) will be saved. This needs to be done
    // BEFORE other methods write arrays to HDF5.

    initHDF5Groups();

    // Allocate memory for the different maps

    subPixelMap.zeros(numRowsSubPixelMap, numColumnsSubPixelMap);
    flatfieldMap.ones((numsubsubfieldsx * (numRowsPixelMap - 2*overlapx) + 2*overlapx) * numSubPixelsPerPixel, (numsubsubfieldsy * (numColumnsPixelMap - 2*overlapy) + 2*overlapy) * numSubPixelsPerPixel);

    if (includeFlatfield)
    {
        // Generate the flatfield map

        generateFlatfieldMap();
    }

    if (includeBFE)
    {
        // Generate Guyonnet coefficients

        generateGuyonnetCoefficients();
    }

    // Initialize and load the PSF. This will open the PSF HDF5 file and perform some basic checking,
    // Then select the proper PSF for the given subfield. Should only be done after calling configure().

    psf = new SymmetricalPointSpreadFunction(configParam, hdf5file);
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
DetectorWithSymmetricalMappedPSF::~DetectorWithSymmetricalMappedPSF()
{
    flushOutput();
    delete psf;
}

/**
 * \brief Configure the DetectorWithSymmetricalMappedPSF object using the given
 *        configuration parameters.
 * 
 * \param configParam: Configuration parameters.
 **/
void DetectorWithSymmetricalMappedPSF::configure(ConfigurationParameters &configParam)
{
    flatfieldNoiseRMS = configParam.getDouble("CCD/FlatfieldNoiseRMS");
    includeFlatfield = configParam.getBoolean("CCD/IncludeFlatfield");
    includeConvolution = configParam.getBoolean("CCD/IncludeConvolution");

    writeSubPixelImagesToHDF5 = configParam.getBoolean(
        "ControlHDF5Content/WriteSubPixelImages");

    numSubPixelsPerPixel = configParam.getInteger("SubField/SubPixels");

    // Configuration parameters for the noise source random seeds

    flatfieldSeed = configParam.getLong("RandomSeeds/FlatFieldSeed");

    // Treat the specific configurations for a Mapped PSF

    string psfModel = configParam.getString("PSF/Model");

    if ((psfModel == "MappedGaussian") || (psfModel == "MappedFromFileSymmetrical"))
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
}

/**
 * \brief Set the PSF map for the sub-field.
 */
void DetectorWithSymmetricalMappedPSF::setPsfForSubfield(int subsubfieldx, int subsubfieldy)
{
    // There is one PSF for the entire subfield, which we take the one of the center
    // of the subfield.

    double xFPmm, yFPmm;
    tie(xFPmm, yFPmm) = getFocalPlaneCoordinatesOfSubfieldCenter(subsubfieldx, subsubfieldy);
    int fieldnumber = subsubfieldx * numsubsubfieldsy + subsubfieldy; //%% unique number for subsubfield
    int fieldmax = numsubsubfieldsx * numsubsubfieldsy -1; //%% total number of subsubfields -1 for counting start at 0

    // Get the 'user specified' angular distance to the optical axis from the psf.
    // If the user didn't specify an angular distance, calculate it from the given
    // focal plane coordinates.

    double radius = psf->getRequestedDistanceToOpticalAxis();

    if (radius < 0.0)
    {
        radius = camera.getGnomonicRadialDistanceFromOpticalAxis(xFPmm, yFPmm);
    }

    psf->select(radius, fieldnumber, fieldmax);

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
    psf->rotate(angle, fieldnumber, fieldmax);

    // Rebin the psfMap to the number of sub-pixels per pixel used for the Detector

    psfVector = psf->rebinToSubPixels(numSubPixelsPerPixel, fieldnumber);

    // Allow the convolver to precompute some stuff given the PSF, so that it doesn't
    // need to be recomputed every convolution.
    for (int binnumber=0; binnumber<wave_bins; binnumber++)
    {
        convolver.initialise(numRowsSubPixelMap, numColumnsSubPixelMap, psfVector[binnumber], binnumber, wave_bins, fieldnumber, fieldmax);
    }
}
