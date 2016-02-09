#include "PointSpreadFunction.h"
#include "Units.h"
#include "Exceptions.h"
#include "Logger.h"


PointSpreadFunction::PointSpreadFunction(ConfigurationParameters &configParam)
{
    configure(configParam);

    groupName = "6000";  // this is currently the only group defined in the HDF5 file

    isSelected = false;
    
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








void PointSpreadFunction::select()
{

}







void PointSpreadFunction::rotate(double angle)
{
}





void PointSpreadFunction::rebin()
{

}

