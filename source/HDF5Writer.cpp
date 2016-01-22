
#include "HDF5Writer.h"



HDF5Writer::HDF5Writer(HDF5File &hdf5File)
: hdf5File(hdf5File)
{

}





HDF5Writer::~HDF5Writer()
{

}





// HDF5Writer::initHDF5Groups()
// 
// PURPOSE: Each HDF5Writer class has its own groups in the HDF5 file where it writes
//          information. These groups have to be created once, at the very beginning.
//
// DEFAULT: do nothing

void HDF5Writer::initHDF5Groups()
{

}







// HDF5Writer::flushOutput()
// 
// PURPOSE: before closing the HDF5 file, flushOutput will be called so that
//          all classes deriving from hdf5writer will have the chance to write
//          any remaining information to the HDF5 file.
//
// DEFAULT: do nothing

void HDF5Writer::flushOutput()
{

}
