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
 * \brief Configure the PointSpreadFunction object using the ConfigurationParameters
 * 
 * \param[in] configParam: the configuration parameters 
 **/

void PointSpreadFunction::configure(ConfigurationParameters &cp)
{
    absolutePath = cp.getAbsoluteFilename("PSF/Filename");
    numberOfSubPixelsPerPixel = cp.getInteger("PSF/SubPixels");
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
        psfMap = rotateArray(psfMap, newAngle);

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








