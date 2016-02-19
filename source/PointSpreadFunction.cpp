#include "PointSpreadFunction.h"
#include "Units.h"
#include "Exceptions.h"
#include "Logger.h"
#include "ArrayOperations.h"

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






PointSpreadFunction::~PointSpreadFunction()
{
    hdf5file->close();
    delete hdf5file;
}







void PointSpreadFunction::configure(ConfigurationParameters &cp)
{
    location = cp.getAbsoluteFilename("PSF/Filename");
}







void PointSpreadFunction::select(double radius)
{
    if (isSelected)
    {
        Log.warning("Another PSF was previously selected.");
    }

    // Convert radius into the string id that identifies the psf dataset in the HDF5 file
    // We work with a lookup table psfdata::radius which contains fixed radius values for which PSF data
    // was generated. 
    
    arma::vec rads = psfdata::radius - radius;
    rads = abs(rads);

    arma::uword index;
    double radius_lut = rads.min(index);

    if (index > psfdata::radius.n_elem-1)
    {
        Log.warning("Radius index out of bounds.");
        index = psfdata::radius.n_elem-1;
    }

    string id = "ar" + to_string(int(psfdata::radius(index) * 1000));

    Log.debug("Identifier for selected PSF is " + id);

    // Load the psf array into the psfMap
    
    hdf5file->readArray("/" + groupName, "ar00000", psfMap);

    isSelected = true;
}







void PointSpreadFunction::rotate(double angle)
{
    if (isRotated)
    {
        Log.warning("The PSF has been previously rotated and will not be rotated again because of inaccuracies.");
        Log.warning("TODO: reload the PSF from the inputfile before rotating.");
    }
    else
    {
        psfMap = rotateArray(psfMap, angle);
        isRotated = true;        
    }
}





void PointSpreadFunction::rebin()
{

}

