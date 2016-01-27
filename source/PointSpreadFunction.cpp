#include "PointSpreadFunction.h"
#include "Exceptions.h"
#include "Logger.h"


PointSpreadFunction::PointSpreadFunction(ConfigurationParameters &cp)
{
    loadConfiguration(cp);
}






PointSpreadFunction::~PointSpreadFunction()
{
    hdf5file->close();
    delete hdf5file;
}


void PointSpreadFunction::loadConfiguration(ConfigurationParameters &cp)
{
    location = cp.getAbsoluteFilename("PSF/Filename");
    groupName = "6000";  // this is currently the only group defined in the HDF5 file
    
    hdf5file = new HDF5File(location);

    if ( !hdf5file->hasGroup(groupName) )
    {
        throw FileException("The HDF5 file doesn't contain the expected group \"" + groupName + "\".");
    }

}




void PointSpreadFunction::select()
{

}





void PointSpreadFunction::rotate()
{

}





void PointSpreadFunction::rebin()
{

}

