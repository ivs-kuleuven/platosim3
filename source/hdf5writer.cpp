
#include "hdf5writer.h"



Hdf5Writer::Hdf5Writer(Hdf5File &hdf5File, string groupName)
hdf5File(hdf5File), hdf5GroupName(groupName)
{

}





~Hdf5Writer::Hdf5Writer()
{

}




// Default: do nothing

void flushOutput()
{

}
