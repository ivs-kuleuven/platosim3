/**
 * \class     PointSpreadFunction 
 * 
 * \brief     The PointSpreadFunction provides the PSF at different positions in the field.
 * 
 * \details
 * 
 * The PSF intensity maps are provided as a HDF5 file which contains several groups of PSFs for 
 * different temperature and different positions in the field. Currently only PSFs for 6000K are
 * available. 
 * 
 * 
 * \todo
 * 
 * There are two hardcoded values used in this class, i.e. the top-level groupName of the HDF5 PSF file, 
 * and the name of the dataset containing the PSF array. The latter should actually be generated from the 
 * ra, dec from the center of the sub-field, but this information is currently not passed into the select()
 * method.
 * 
 */

#include "PointSpreadFunction.h"








/**
 * \brief      Constructor
 *
 * \details
 * 
 * Load the configuration parameters and initialize the internal variables. The default group 
 * for which PSF will be loaded is currently set to "6000" (which is the only group available at the moment).
 * 
 * The HDF5 PSF file is loaded and some basic checking is done.
 * 
 * \param      configParam  configuration parameters for the PSF
 */

PointSpreadFunction::PointSpreadFunction(ConfigurationParameters &configParam)
{
    configure(configParam);

    isSelected = false;
    isRotated = false;
    
    if ( ! FileUtilities::fileExists(absolutePath) )
    {
        throw FileException("PointSpreadFunction: trying to load the PSF HDF5 file (" + absolutePath + "), but file doesn't exist.");
    }

    try
    {
        hdf5file.open(absolutePath);
    }
    catch(H5::FileIException ex)
    {
        Log.error("H5::FileIException: " + string(ex.getCDetailMsg()));
        throw H5FileException("PointSpreadFunction: Could not open HDF5 file: " + absolutePath);
    }

    string groupName = "T6000";
    if ( ! hdf5file.hasGroup(groupName) )
    {
        throw H5FileException("PointSpreadFunction: The HDF5 file (" + absolutePath + ") doesn't contain the expected group \"" + groupName + "\".");
    }

    Log.info("PointSpreadFunction: Opened the HDF5 file " + absolutePath + " containing the PSFs");

}






/**
 * \brief      Destructor
 * 
 * \Details
 * 
 * Close the HDF5 file and release the memory.
 * 
 */
PointSpreadFunction::~PointSpreadFunction()
{
    hdf5file.close();
}








/**
 * @brief      Return a two-dimensional Gaussian function
 * 
 * @details
 * 
 * The Gaussian PSF is calculated from the following equation:
 * 
 * \f[
 *   f(x, y)  =  \frac{1}{2\pi\sigma^{2}} e^{-[(x - \mu_x)^{2} + (y - \mu_y)^{2}] / (2\sigma^{2})}
 * \f]
 *
 *
 * @return     a 2D Gaussian PSF
 */
arma::fmat PointSpreadFunction::getGaussianPsf()
{

    double centerRow = numberOfSubPixelsPerPixel * numberOfPixels / 2.0;
    double centerColumn = numberOfSubPixelsPerPixel * numberOfPixels / 2.0;

    if (sigma < 1.0 / numberOfSubPixelsPerPixel)
    {
        throw IllegalArgumentException("The width of the Gaussian PSF must be larger than the size of a subpixel.");
    }

    // Generate the Gaussian PSF at the sub-pixel level

    double width = sigma * numberOfSubPixelsPerPixel;
    double normalizationFactor = 1.0 / (width * width * 2.0 * Constants::PI);
    double denominator = 2.0 * width * width;

    arma::fmat gaussianPsf (numberOfPixels * numberOfSubPixelsPerPixel, numberOfPixels * numberOfSubPixelsPerPixel, arma::fill::zeros);

    for (unsigned int xi = 0; xi < gaussianPsf.n_cols; xi++)
    {
        for (unsigned int xj = 0; xj < gaussianPsf.n_rows; xj++)
        {
            // FIXME: This a equation can probably be optimized with Armadillo functionality
            gaussianPsf(xi, xj) = normalizationFactor 
                * exp( - (pow(xi - centerColumn, 2.0) + pow(xj - centerRow, 2.0)) / denominator);
        }
    }

    return gaussianPsf;
}


/**
 * \brief Configure the PointSpreadFunction object using the ConfigurationParameters
 * 
 * \param[in] configParam: the configuration parameters 
 **/

void PointSpreadFunction::configure(ConfigurationParameters &cp)
{
    isGaussian                = cp.getBoolean("PSF/UseGauss");
    absolutePath              = cp.getAbsoluteFilename("PSF/Filename");
    numberOfSubPixelsPerPixel = cp.getInteger("PSF/NumberOfSubPixels");
    numberOfPixels            = cp.getInteger("PSF/NumberOfPixels");
    sigma                     = cp.getDouble("PSF/Sigma");
}







/**
 * \brief      Select the proper PSF matching the given radius closest.
 *
 * \details
 * 
 * The generated PSF is position dependent and a look-up table has been provided that contains the 
 * field radial coordinate along the line of sight for which the PSF was generated.
 *  
 * The appropriate PSF will be selected, i.e. the PSF for which the angular distance of the centre 
 * of the sub-field from the centre of the focal plane matches best.
 * 
 * \param[in]  radius  angular separation of the source for which to select the PSF [radians]
 */
void PointSpreadFunction::select(double radius)
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

    if ( ! hdf5file.hasGroup(groupName) )
    {
        throw FileException("PointSpreadFunction: The HDF5 file (" + absolutePath + ") doesn't contain the expected group \"" + groupName + "\".");
    }

    if ( ! hdf5file.hasDataset(groupName, azimuthDataset) )
    {
        throw FileException("PointSpreadFunction: The HDF5 file (" + absolutePath + ") doesn't contain the expected dataset \"" + azimuthDataset + "\".");
    }

    // Load the psf array into the psfMap
    
    hdf5file.readArray("/" + groupName, azimuthDataset, psfMap);
    
    // The PSFs that are currently used are rotated with respect to the focal plane x-axis.
    // The rotation angle is given as an attribute to the dataset that contains the PSF.

    double angle = hdf5file.readAttribute(groupName, azimuthDataset, "orientation");

    rotationAngle = deg2rad(angle);

    Log.debug("PointSpreadFunction: Selected PSF " + groupName + "/" + azimuthDataset + ", rotation set to " + dtos(angle) + " degrees.");

    isSelected = true;
}







/**
 * \brief      Rotate the PSF with the given angle.
 * 
 * \details
 * 
 * Beware that the PSF that has been provided is already rotated, this will be taken into account.
 *
 * \param[in]  angle  angle by which the PSF should be rotated [radians]
 */
void PointSpreadFunction::rotate(double angle)
{

    // We do not need to rotate a Gaussian PSF

    if (isGaussian)
        return;

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

        Log.debug("PointSpreadFunction: rotated current PSF over angle " + to_string(rad2deg(newAngle)) + " deg");    
    }
}







/**
 * @brief      return a reference to the Armadillo array that contains the selected PSF
 *
 * @return     a float array with the PSF
 */
arma::Mat<float> PointSpreadFunction::getPsfMap()
{
    return psfMap;
}






/**
 * @brief      Rebin the PSF map to the target number of subpixels
 *
 * @details
 * 
 * The number of subpixels used to generate the PSF is not necessarily the same
 * as the number of subpixels per pixel for the detector. So, the PSF needs to be
 * rebinned to the number of subpixels per pixel for the detector, which is given
 * as targetSubPixels.
 * 
 * @param[in]  targetSubPixels  the target number of subpixels
 */
void PointSpreadFunction::rebin(unsigned int targetSubPixels)
{
    unsigned int binSize = psfMap.n_rows / numberOfSubPixelsPerPixel * targetSubPixels;

    psfMap = ArrayOperations::rebin(psfMap, binSize, binSize);

    isRebinned = true;
}




