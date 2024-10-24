
// Code that can be used to test HDF5 stuff
//
// Compile with
// clang++ -o testhdf5 testhdf5.cpp -I../../dependencies/Installs/hdf5-1.10.2/include  -L../../dependencies/Installs/hdf5-1.10.2/lib -lhdf5 -lhdf5_cpp -stdlib=libc++ -std=c++14
//

#include <iostream>
#include <string>
#include "H5Cpp.h"

using namespace std;

void testFunc1()
{
    H5::Group group;
    H5::DataSet dataset;
    H5::Attribute attr;

    string fileName = "/Users/joris/Downloads/blueRealPSF.hdf5";
    H5::H5File *file = new H5::H5File(fileName.c_str(), H5F_ACC_RDONLY);

    string groupName = "/";
    try
    {
        group = file->openGroup(groupName.c_str());
    }
    catch (H5::Exception ex) {}

    string datasetName = "1";
    try
    {
        dataset = group.openDataSet(datasetName.c_str());
    }
    catch (H5::Exception ex) {}

    string attributeName = "starPointing1";
    try 
    {  
        attr = dataset.openAttribute(attributeName.c_str());
    }
    catch (H5::AttributeIException error) {}

    double value = 0.0;

    H5::DataType type = attr.getDataType();
    attr.read(type, &value);

    cout << "Attribute value: " << value << endl;

    file->close();
    delete file;
}





void testFunc2()
{
    H5::Group group;
    H5::DataSet dataset;
    H5::Attribute attr;

    string fileName = "/Users/joris/Downloads/blueRealPSF.hdf5";
    H5::H5File *file = new H5::H5File(fileName.c_str(), H5F_ACC_RDONLY);

    string groupName = "/";
    try
    {
        group = file->openGroup(groupName.c_str());
    }
    catch (H5::Exception ex) {}

    string datasetName = "1";
    try
    {
        dataset = group.openDataSet(datasetName.c_str());
    }
    catch (H5::Exception ex) {}

    string attributeName = "centerCoordinates1";
    try 
    {  
        attr = dataset.openAttribute(attributeName.c_str());
    }
    catch (H5::AttributeIException error) { cout << "Attribute exception" << endl;}


    H5::StrType stringType = attr.getStrType();

    string attributeValue;
    attr.read(stringType, attributeValue);

    cout << "Attribute value: |" << attributeValue << "|" << endl;

    double myDouble = stod(attributeValue); 
    cout << "myDouble = " << myDouble << endl;

    file->close();
    delete file;
}













int main()
{
    testFunc1();
    testFunc2();
    return 0;
}
