
#include "hdf5writer.h"



HDF5Writer::HDF5Writer(HDF5File &hdf5File, string groupName)
hdf5File(hdf5File), hdf5GroupName(groupName)
{

}





~HDF5Writer::HDF5Writer()
{

}




// Default: do nothing

void flushOutput()
{

}
