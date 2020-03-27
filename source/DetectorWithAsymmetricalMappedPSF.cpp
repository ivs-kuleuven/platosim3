#include "DetectorWithAsymmetricalMappedPSF.h"

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

DetectorWithAsymmetricalMappedPSF::DetectorWithAsymmetricalMappedPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
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
    flatfieldMap.ones(numRowsSubPixelMap, numColumnsSubPixelMap);

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

    psf = new AsymmetricalPointSpreadFunction(configParam, hdf5file);
    setPsfForSubfield();
}









/**
 * Destructor.
 */
DetectorWithAsymmetricalMappedPSF::~DetectorWithAsymmetricalMappedPSF()
{
    flushOutput();
    delete psf;
}








/**
 * \brief Configure the DetectorWithAsymmetricalMappedPSF object using the given
 *        configuration parameters.
 * 
 * \param configParam: Configuration parameters.
 */
void DetectorWithAsymmetricalMappedPSF::configure(ConfigurationParameters &configParam)
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

    if (psfModel == "MappedFromFileAsymmetrical")
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
 * 
 * \details The PSF that is selected is dependent on the user input.
 */
void DetectorWithAsymmetricalMappedPSF::setPsfForSubfield()
{
    // There is one PSF for the entire subfield, which we take the one of the center
    // of the subfield.

    double xFPmm, yFPmm;
    tie(xFPmm, yFPmm) = getFocalPlaneCoordinatesOfSubfieldCenter();

    psf->select(xFPmm, yFPmm);

    //  Compensate for the orientation of the CCD wrt focal plane orientation.

    psf->rotate(-orientationAngle);

    // Rebin the psfMap to the number of sub-pixels per pixel used for the Detector

    psfMap = psf->rebinToSubPixels(numSubPixelsPerPixel);

    // Allow the convolver to precompute some stuff given the PSF, so that it doesn't
    // need to be recomputed every convolution.

    convolver.initialise(numRowsSubPixelMap, numColumnsSubPixelMap, psfMap);
}
