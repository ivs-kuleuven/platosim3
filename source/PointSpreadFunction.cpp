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
#include "Units.h"
#include "Exceptions.h"
#include "Logger.h"
#include "ArrayOperations.h"









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
    
    hdf5file = new HDF5File(location);

    if ( !hdf5file->hasGroup("T6000") )
    {
        throw FileException("The HDF5 file doesn't contain the expected group \"" + groupName + "\".");
    }

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
    hdf5file->close();
    delete hdf5file;
}






/**
 * \brief Configure the PointSpreadFunction object using the ConfigurationParameters
 * 
 * \param[in] configParam: the configuration parameters 
 **/

void PointSpreadFunction::configure(ConfigurationParameters &cp)
{
    location = cp.getAbsoluteFilename("PSF/Filename");
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
    // TODO: Should we take any action if different PSFs are selected for this object?

    if (isSelected)
    {
        Log.warning("Another PSF was previously selected.");
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
        Log.warning("Radius index out of bounds.");
        index = psfdata::radius.n_elem-1;
    }

    string angularRadiusGroup = "ar" + to_string(int(psfdata::radius(index) * 1000));
    string temperatureGroup = "T6000";  // TODO: hardcoded value! 
    string azimuthDataset = "az0";      // TODO: hardcoded value!

    string groupName = temperatureGroup + "/" + angularRadiusGroup;

    if ( ! hdf5file->hasGroup(groupName) )
    {
        throw FileException("The HDF5 file doesn't contain the expected group \"" + groupName + "\".");
    }

    if ( ! hdf5file->hasDataset(groupName, azimuthDataset) )
    {
        throw FileException("The HDF5 file doesn't contain the expected dataset \"" + azimuthDataset + "\".");
    }

    // Load the psf array into the psfMap
    
    hdf5file->readArray("/" + groupName, azimuthDataset, psfMap);
    
    // The PSFs that are currently used are rotated with respect to the focal plane x-axis.
    // The rotation angle is given as an attribute to the dataset that contains the PSF.

    double angle = hdf5file->readAttribute(groupName, azimuthDataset, "orientation");
    Log.debug("PointSpreadFunction::select: orientation = " + to_string(angle));

    rotationAngle = deg2rad(angle);

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
        Log.warning("The PSF has been previously rotated and will not be rotated again because of inaccuracies.");
        Log.warning("TODO: reload the PSF from the inputfile before rotating.");
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
    }
}





void PointSpreadFunction::rebin()
{

}

