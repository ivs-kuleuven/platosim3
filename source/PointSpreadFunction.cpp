#include "PointSpreadFunction.h"

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




bool hasDataset(const string &)
{
    return false;
}


void PointSpreadFunction::loadConfiguration(ConfigurationParameters &cp)
{
    location = cp.getAbsoluteFileName("PSF/FileName");
    groupName = "6000";  // this is currently the only group defined in the HDF5 file
    
    hdf5file = new HDF5File(location);

    if ( !hdf5file->hasGroup(groupName) )
    {
        Log.error("PointSpreadFunction: The HDF5 file doesn't contain the expected group \"" + groupName + "\".");
        exit(1);
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

