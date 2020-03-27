/**
 * \class SymmetricalPointSpreadFunction 
 * 
 * \brief Class for rotationally symmetrical PSFs, at different location in the field (i.e. at
 *        different angular distances from the optical axis).
 * 
 * \details The PSF intensity maps are provided as an HDF5 file, which contains one group per temperature, with
 *          sub-groups for different positions in the field (i.e. different angular distances from the optical axis).
 *          Currently only PSFs for a temperature of 6000K are available.
 * 
 * 
 * \todo There are two hardcoded values used in this class, i.e. the top-level group name of the HDF5 PSF file
 *       and the name of the dataset containing the PSF array. The latter should actually be generated from the 
 *       ra, dec from the center of the sub-field, but this information is currently not passed into the select()
 *       method.
 */

#include "SymmetricalPointSpreadFunction.h"





/**
 * \brief Constructor.
 *
 * \details Loads the configuration parameters and initialises the internal variables. The default group 
 *          for which PSF will be loaded is currently set to "T6000" (which is the only group available
 *          at the moment).  The HDF5 PSF file is loaded and some basic checking is done.
 * 
 * \param configParam: Configuration parameters for the PSF.
 * 
 * \param hdf5file: HDF5 file from which to read the PSF.
 */
SymmetricalPointSpreadFunction::SymmetricalPointSpreadFunction(ConfigurationParameters &configParam, HDF5File &hdf5file) : PointSpreadFunction(configParam, hdf5file)
{
    // Create the groups in the HDF5 file where the different PSFs and their description will be saved.

    initHDF5Groups();

    // Parse the parameters from the configuration file.

    configure(configParam);

    isSelected = false;
    isRotated = false;
    
    if (isLoadedFromFile)
    {
        // Prepare the psfFile by performing some basic checks

        if (!FileUtilities::fileExists(absolutePath))
        {
            throw FileException("SymmetricalPointSpreadFunction: trying to load the PSF HDF5 file (" + absolutePath + "), but file doesn't exist.");
        }
    
        try
        {
            psfFile.open(absolutePath);
        }
        catch(H5::FileIException ex)
        {
            Log.error("H5::FileIException: " + string(ex.getCDetailMsg()));
            throw H5FileException("SymmetricalPointSpreadFunction: Could not open HDF5 file: " + absolutePath);
        }
    
        string groupName = "T6000";
        if ( ! psfFile.hasGroup(groupName) )
        {
            throw H5FileException("SymmetricalPointSpreadFunction: The HDF5 file (" + absolutePath + ") doesn't contain the expected group \"" + groupName + "\".");
        }
    
        Log.info("SymmetricalPointSpreadFunction: Opened the HDF5 file " + absolutePath + " containing the PSFs");
    }
}





/**
 * \brief Returns a two-dimensional Gaussian function.
 *
 * \details The Gaussian PSF is calculated from the following equation:
 *
 *             \f[
 *                f(x, y)  =  \frac{1}{2\pi\sigma^{2}} e^{-[(x - \mu_x)^{2} + (y - \mu_y)^{2}] / (2\sigma^{2})}
 *             \f]
 *            
 *             where \f$\mu_x\f$ and \f$\mu_y\f$ are the mean (center points) and 
 *             \f$\sigma\f$ is the standard deviation.
 *            
 *             The Gaussian PSF is centered in an array at subpixel level, i.e. the array
 *             has dimension [numberOfSubPixelsPerPixel * numberOfPixels,
 *             numberOfSubPixelsPerPixel * numberOfPixels], but the standard deviation 
 *             \f$\sigma\f$ is defined at pixel level.
 *
 * \return: 2D Gaussian PSF.
 */
arma::fmat SymmetricalPointSpreadFunction::getGaussianPsf()
{
    double centerRow = numberOfSubPixelsPerPixel * numberOfPixels / 2.0;
    double centerColumn = numberOfSubPixelsPerPixel * numberOfPixels / 2.0;

    if (sigma < 1.0 / numberOfSubPixelsPerPixel)
    {
        throw IllegalArgumentException("The width of the Gaussian PSF must be larger than the size of a subpixel.");
    }

    // Generate the Gaussian PSF at the sub-pixel level

    double width = sigma * numberOfSubPixelsPerPixel;
    double denominator = 2.0 * width * width;
    double sumPSF = 0.0;

    arma::fmat gaussianPsf (numberOfPixels * numberOfSubPixelsPerPixel, numberOfPixels * numberOfSubPixelsPerPixel, arma::fill::zeros);

    for (unsigned int xi = 0; xi < gaussianPsf.n_cols; xi++)
    {
        for (unsigned int xj = 0; xj < gaussianPsf.n_rows; xj++)
        {
            gaussianPsf(xj, xi) = exp( - (pow(xi - centerColumn, 2) + pow(xj - centerRow, 2)) / denominator);
            sumPSF += gaussianPsf(xj, xi);
        }
    }

    // Normalize the gaussian, so that the flux is conserved
    
    gaussianPsf /= sumPSF;

    return gaussianPsf;
}





/**
 * \brief Configures the SymmetricalPointSpreadFunction object using the given
 *        configuration parameters
 *
 * \param[in] configParam: Configuration parameters.
 */
void SymmetricalPointSpreadFunction::configure(ConfigurationParameters &configParam)
{
   string model = configParam.getString("PSF/Model");

    if (model == "MappedGaussian")
    {
        // The user specified to use a Gaussian shape PSF
        // The number of sub-pixels per pixel that will be used to calculate the 
        // Gaussian PSF is equal to the number of sub-pixels per pixels for the sub-field.

        isGaussian  = true;
        sigma = configParam.getDouble("PSF/MappedGaussian/Sigma");
        numberOfPixels = configParam.getInteger("PSF/MappedGaussian/NumberOfPixels");

        // The Gaussian PSF shall be created with a resolution equal to that of the sub-field
        
        numberOfSubPixelsPerPixel = configParam.getInteger("SubField/SubPixels");
    } 
    else if (model == "MappedFromFileSymmetrical")
    {
        // The user specified to use the pre-calculated PSFs from file
        // The number of sub-pixels per pixel is derived from the number of pixels specified 
        // and the size of the array in the file. (The number of sub-pixels should be in the file).

        isLoadedFromFile = true;
        absolutePath = configParam.getAbsoluteFilename("PSF/MappedFromFileSymmetrical/Filename");
        numberOfPixels = configParam.getInteger("PSF/MappedFromFileSymmetrical/NumberOfPixels");
        requestedDistanceToOA = deg2rad(configParam.getDouble("PSF/MappedFromFileSymmetrical/DistanceToOA"));
        requestedRotationAngle = deg2rad(configParam.getDouble("PSF/MappedFromFileSymmetrical/RotationAngle"));
    }
    else
    {
        string errorMessage = "SymmetricalPointSpreadFunction: Model '" + model + "' is not supported.";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
    }
}





/**
 * \brief Returns the distance to the optical axis as requested by the user.  When the returned value 
 *        is negative [-1], the user input will be ignored
 *
 * \return Distance to optical axis [rad].
 */
double SymmetricalPointSpreadFunction::getRequestedDistanceToOpticalAxis()
{
    return requestedDistanceToOA;
}





/**
 * \brief Returns the orientation angle of the PSF as requested by the user. This orientation is
 *        what the user specified in the input file and is not necessarily equal to the rotation 
 *        angle of the PSF (which is also compensated for the CCD orientation).
 *
 * \return Orientation as requested by the user.
 */
double SymmetricalPointSpreadFunction::getRequestedRotationAngle()
{
    return requestedRotationAngle;
}





/**
 * \brief Selects the proper PSF matching the given radius closest.  The generated PSF is position dependent 
 *        and a look-up table has been provided that contains the field radial coordinate along the line-of-sight
 *        for which the PSF was generated.  The appropriate PSF will be selected, i.e. the PSF for which the angular 
 *         distance of the centre of the sub-field from the centre of the focal plane matches best.
 * 
 * \param[in] radius: Angular distance of the source for which to select the PSF from the optical axis [radians].
 */
void SymmetricalPointSpreadFunction::select(double radius)
{
    using StringUtilities::dtos;
    
    // TODO: Should we take any action if different PSFs are selected for this object?

    if (isSelected)
    {
        Log.warning("PointSpreadFunction: Another PSF was previously selected.");
    }

    if (isGaussian)
    {
        psfMap = getGaussianPsf();
        isSelected = true;
        rotationAngle = 0.0;

        hdf5File.writeAttribute("/PSF", "selectedPSF", "Gaussian PSF selected with sigma=" + to_string(sigma));

        return;
    }

    radius = rad2deg(radius);

    // Convert radius into the string angularRadiusGroup that identifies the psf dataset in the HDF5 file.
    // We work with a lookup table psfdata::radius which contains fixed radius values for which PSF data
    // was generated. The algorithm is to select the PSF with radius closest to the given radius by 
    // subtracting the given radius from the tabulated radius and then selecting the lowest value to find the index.
    
    arma::vec rads = psfdata::radius - radius;
    rads = abs(rads);

    arma::uword index;
    rads.min(index);

    if (index > psfdata::radius.n_elem-1)
    {
        Log.warning("PointSpreadFunction: Radius index (" + to_string(index) + ") is out of bounds.");
        index = psfdata::radius.n_elem-1;
    }

    stringstream myStream;
    myStream << "ar" << setfill('0') << setw(5) << int(psfdata::radius(index) * 1000);

    const string angularRadiusGroup = myStream.str();
    const string temperatureGroup = "T6000";  // TODO: hardcoded value! 
    const string azimuthDataset = "az0";      // TODO: hardcoded value!

    string groupName = temperatureGroup + "/" + angularRadiusGroup;

    if ( ! psfFile.hasGroup(groupName) )
    {
        throw FileException("PointSpreadFunction: The HDF5 file (" + absolutePath + ") doesn't contain the expected group \"" + groupName + "\".");
    }

    if ( ! psfFile.hasDataset(groupName, azimuthDataset) )
    {
        throw FileException("PointSpreadFunction: The HDF5 file (" + absolutePath + ") doesn't contain the expected dataset \"" + azimuthDataset + "\".");
    }

    // Load the psf array into the psfMap
    
    psfFile.readArray("/" + groupName, azimuthDataset, psfMap);
    
    // The PSFs that are currently used are rotated with respect to the focal plane x-axis.
    // The rotation angle is given as an attribute to the dataset that contains the PSF.

    double angle = psfFile.readDoubleDatasetAttribute(groupName, azimuthDataset, "orientation");

    rotationAngle = deg2rad(angle);

    // We should be able to read the number of sub-pixels per pixel that was used to generate the PSFs
    // from an attribute in the HDF5 file. Unfortunately, this is not available and we therefore derive 
    // the number from the array size of the psfMap and the number of pixels, currently specified 
    // in the input file.

    numberOfSubPixelsPerPixel = psfMap.n_rows / numberOfPixels;

    Log.debug("PointSpreadFunction: Selected PSF " + groupName + "/" + azimuthDataset + ", rotation set to " + dtos(angle) + " degrees.");

    hdf5File.writeAttribute("/PSF", "selectedPSF", "Realistic PSF selected from group " + groupName + "/" + azimuthDataset + ".");

    isSelected = true;
}





/**
 * \brief Rotates the PSF with the given angle.  Beware that the PSF that has been provided is already 
 *        rotated, this will be taken into account.
 *
 * \param[in] angle: Angle by which the PSF should be rotated [radians].
 */
void SymmetricalPointSpreadFunction::rotate(double angle)
{
    // We do not need to rotate a Gaussian PSF
    // Even if we do not rotate the Gaussian PSF, we do save the psfMap as a rotatedPSF.
    // This is to keep consistency in the output file where we do not save the selected PSF,
    // but we do save the rotated PSF.

    if (isGaussian)
    {
        hdf5File.writeArray("/PSF", "rotatedPSF", psfMap);
        hdf5File.writeAttribute("/PSF", "rotationAngle", rotationAngle);
        return;
    }

    if (isRotated)
    {
        Log.warning("PointSpreadFunction: Ignoring rotation: PSF was already rotated before.");
    }
    else
    {
        // This part of the code currently assumes that the PSF has not been rotated yet.
        // The rotationAngle of the generated PSFs is 45 degrees, so the actual rotation should be
        // the requested angle minus the rotationAngle of the PSF.

        double newAngle = angle - rotationAngle;
        psfMap = ArrayOperations::rotateArray(psfMap, newAngle);

        rotationAngle = newAngle;
        isRotated = true;    

        psfMap /= arma::accu(psfMap);

        Log.debug("PointSpreadFunction: rotated current PSF over angle " + to_string(rad2deg(newAngle)) + " deg");

        // Write the psfMap of the rotated PSF to the HDF5 output file

        hdf5File.writeArray("/PSF", "rotatedPSF", psfMap);
        hdf5File.writeAttribute("/PSF", "rotationAngle", rotationAngle);
    }
}
