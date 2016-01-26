#include "PointSpreadFunction.h"


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
    location = cp.getAbsoluteFileName("Camera/PSFFileName");
    groupName = "6000";  // this is currently the only group defined in the HDF5 file
    
    hdf5file = new HDF5File(location);

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

