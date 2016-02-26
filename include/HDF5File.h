
#ifndef HDF5FILE_H
#define HDF5FILE_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "H5Cpp.h"
#include "memory.h"
#include "armadillo"

#include "Logger.h"
#include "HDF5Exceptions.h"

using namespace std;



class HDF5File
{
    public:

        HDF5File();
        HDF5File(string filename, bool readonly=true);
        ~HDF5File();

        void open(string filename, bool readonly=true);
        void close();
        
        bool hasGroup(string groupName);
        void createGroup(string groupName);

        bool hasDataset(string groupName, string datasetName);

        void writeAttribute(string groupName, string attributeName, string attributeValue);
        void writeAttribute(string groupName, string attributeName, int attributeValue);
        void writeAttribute(string groupName, string attributeName, long attributeValue);
        void writeAttribute(string groupName, string attributeName, double attributeValue);

        void writeArray(string groupName, string arrayName, int*    array, int size);
        void writeArray(string groupName, string arrayName, float*  array, int size);
        void writeArray(string groupName, string arrayName, double* array, int size);
        void writeArray(string groupName, string arrayName, arma::Mat<float>& A);

        double readAttribute(string groupName, string attributeName);
        double readAttribute(string groupName, string datasetName, string attributeName);

        void readArray(string groupName, string arrayName, arma::Mat<float>& A);


    protected:


    private:

        H5::H5File *file;
        bool fileIsOpen;

};



bool fileExists(string filename);


#endif
