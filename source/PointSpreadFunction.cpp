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

    groupName = "6000";  // this is currently the only group defined in the HDF5 file

    isSelected = false;
    isRotated = false;
    
    hdf5file = new HDF5File(location);

    if ( !hdf5file->hasGroup(groupName) )
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
    if (isSelected)
    {
        Log.warning("Another PSF was previously selected.");
    }

    radius = rad2deg(radius);
    
    // Convert radius into the string id that identifies the psf dataset in the HDF5 file.
    // We work with a lookup table psfdata::radius which contains fixed radius values for which PSF data
    // was generated. The algorithm is to select the PSF with radius closest to the given radius by 
    // the given radius from the tabulated radius and then selecting the lowest value to find the index.
    
    arma::vec rads = psfdata::radius - radius;
    rads = abs(rads);

    arma::uword index;
    rads.min(index);

    if (index > psfdata::radius.n_elem-1)
    {
        Log.warning("Radius index out of bounds.");
        index = psfdata::radius.n_elem-1;
    }

    string id = "ar" + to_string(int(psfdata::radius(index) * 1000));

    Log.debug("Identifier for selected PSF is " + id);

    // Load the psf array into the psfMap
    
    hdf5file->readArray("/" + groupName, "ar00000", psfMap);
    
    // The PSFs that are currently used are rotated by 45 degrees.
    // So this values is currently hard-coded because it is not provided as part of the delivery pack.
    rotationAngle = 45.0;

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
    // Angles in this method are all in degrees

    angle = rad2deg(angle);

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

