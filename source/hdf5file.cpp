
#include "hdf5file.h"



// Default constructor

HDF5File::HDF5File()
: file(NULL), fileIsOpen(false)
{

}






// Open file constructor
//
// PURPOSE: open an HDF5 file. If it doesn't already exist, create it.
//
// INPUT: - filename: full path of the file
//        - overwrite: true if previously stored information in the file should be erased,
//                     false otherwise

HDF5File::HDF5File(string filename, bool overwrite)
: file(NULL), fileIsOpen(false)
{
    open(filename, overwrite);
}







// Destructor

HDF5File::~HDF5File()
{
    close();
    delete file;
}







// HDF5File::open()
//
// PURPOSE: open an HDF5 file. If it doesn't already exist, create it.
//
// INPUT: - filename: full path of the file
//        - overwrite: true if previously stored information in the file should be erased,
//                     false otherwise
//
// OUTPUT: None

void HDF5File::open(string filename, bool overwrite)
{
    // If there was already a file opened, complain and exit

    if (fileIsOpen)
    {
        Log.error("HDF5File::open(): file " + filename + " already open.");
        close();
        exit(1);
    }

    // Check if the file already exists.
    // If so, use append mode, if not

    if (fileExists(filename))
    {
        if (overwrite)
        {
            // Open file, and erase (truncate) all data previously stored in the file

            file = new H5::H5File(filename.c_str(), H5F_ACC_TRUNC);
        }
        else
        {
            // Open existing file for read/write.
            // Conserve the data previously stored in the file.

            file = new H5::H5File(filename.c_str(), H5F_ACC_RDWR);
        }
    }
    else
    {
        // The file doesn't exist yet. Create a new one.

        file = new H5::H5File(filename.c_str(), H5F_ACC_TRUNC);
    }

    fileIsOpen = true;
    return;
}









// HDF5File::close()
//
// PURPOSE: close the HDF5 file in case it was open. 
//          It's safe to close a file that wasn't even open.
//
// INPUT: None
//
// OUTPUT: None

void HDF5File::close()
{
    if (fileIsOpen)
    {
        file->close();
        delete file;
        file = NULL;
        fileIsOpen = false;
    }
}













// HDF5File::fileExists()
//
// PURPOSE: Check if file exists. This is done by trying to open the file,
//          and verifying if the "good" flag is raised.
//
// INPUT: filename: full path of the file
//
// OUTPUT: true if the file exists, false otherwise

bool HDF5File::fileExists(string filename)
{
    ifstream myFile(filename.c_str()); 
 
    if (myFile.good())
    {
        myFile.close();
        return true;
    } 
    else 
    {
        myFile.close();
        return false;
    }   
}











// HDF5File::createGroup()
//
// PURPOSE: Create a group in an HDF5 file.
//
// INPUT: groupName: full path of the group. Should always start with "/".
//
// OUTPUT: None
// 
// EXAMPLE: createGroup("/my/path/to/subgroup1") will create "subgroup1"
//          in the parent group "/my/path/to". 
//          Note that the parent group is assumed to already exist. 
//          If not, an error will be given. So, to create "/my/path/to/subgroup1" 
//          from scratch, you should call:
//             createGroup("/my");
//             createGroup("/my/path");
//             createGroup("/my/path/to");
//             createGroup("/my/path/to/subgroup1");

void HDF5File::createGroup(string groupName)
{
    // Make sure the path of the group starts with a "/" (i.e. the root folder)

    if (!groupName.compare(0, 1, "/"))
    {
        groupName.insert(0, "/");
    }

    // Find the parent group. 
    // E.g. if the groupName is "/my/path/to/group1" then the parent group is 
    //      "/my/path/to" and the subgroup is "group1".

    auto position = groupName.find_last_of("/");
    string parentGroupName = groupName.substr(0, position);
    string subGroupName = groupName.substr(position+1);

    // If the parent group name is empty, it means that the group should be 
    // placed in the root group. In that case, make sure that the parentGroupName
    // is "/".

    if (parentGroupName.empty())
    {
        parentGroupName = "/";
    }

    // Open the parent group

    H5::Group parentGroup = file->openGroup(parentGroupName.c_str());

    // Create the new subgroup inside the parent group.

    H5::Group subGroup = parentGroup.createGroup(groupName.c_str());
}












// HDF5File::writeAttribute()  for string-valued attributes
//
// PURPOSE: add a string attribute to the HDF5 file.
//          The attribute is stored in "<GroupName>/<attributeName>".
//          
// INPUT: groupName:      string containing the full path of an existing group. Starts with "/".
//        attributeName:  string containing the name of the attribute
//        attributeValue: string containing the attribute value
//
// OUTPUT: None
//
// EXAMPLE: writeAttribute("/InputParameters/ObservingParameters", "CCDPredefinedPosition", "User")
//

void HDF5File::writeAttribute(string groupName, string attributeName, string attributeValue)
{
   // Complain if the file was not first opened
    
    if (!fileIsOpen)
    {
        Log.error("HDF5File::writeAttribute(): file is not open.");
        exit(1);
    }

    // Open the proper group where the input parameter belongs

    H5::Group group = file->openGroup(groupName.c_str());

    // Check whether the attribute is not already in the group. If so, complain.
    // The only way to do so, seems to try to access it and catch the exception
    // if it does not yet exist.

    bool attributeIsAlreadyInGroup = true;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the attribute

        H5::Attribute attr = group.openAttribute(attributeName.c_str());
    }
    catch (H5::AttributeIException error)
    {
        attributeIsAlreadyInGroup = false;
    }

    if (attributeIsAlreadyInGroup)
    {
        Log.error("HDF5File::writeAttribute(): attribute " + attributeName + " already in group " + groupName);
        exit(1);
    }


    // Create and write the attribute to the group
    // H5T_VARIABLE refers to a variable-length string
    // H5S_SCALAR refers to a string scalar attribute

    H5::StrType variableLengthStringType(0, H5T_VARIABLE);
    H5::DataSpace attributeSpace(H5S_SCALAR);
    H5::Attribute attribute = group.createAttribute(attributeName.c_str(), variableLengthStringType, attributeSpace);
    attribute.write(variableLengthStringType, H5std_string(attributeValue.c_str()));

    // That's it

    attribute.close();
    group.close();
    
    return;
}











// HDF5File::writeAttribute() for integer attributes
//
// PURPOSE: add an integer attribute to the HDF5 file.
//          The attribute is stored in "<GroupName>/<attributeName>".
//          
// INPUT: groupName:      string containing the full path of an existing group. Starts with "/".
//        attributeName:  string containing the name of the attribute
//        attributeValue: integer containing the attribute value
//
// OUTPUT: None
//
// EXAMPLE: writeAttribute("/InputParameters/CCD", "CCDSizeX", 4510)
//

void HDF5File::writeAttribute(string groupName, string attributeName, int attributeValue)
{
    // Complain if the file was not first opened
    
    if (!fileIsOpen)
    {
        Log.error("HDF5File::writeAttribute(): file is not open.");
        exit(1);
    }

 
    // Open the proper group where the input parameter belongs

    H5::Group group = file->openGroup(groupName.c_str());

    // Check whether the attribute is not already in the group. If so, complain.
    // The only way to do so, seems to try to access it and catch the exception
    // if it does not yet exist.

    bool attributeIsAlreadyInGroup = true;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the attribute

        H5::Attribute attr = group.openAttribute(attributeName.c_str());
    }
    catch (H5::AttributeIException error)
    {
        attributeIsAlreadyInGroup = false;
    }

    if (attributeIsAlreadyInGroup)
    {
        Log.error("HDF5File::writeAttribute(): attribute " + attributeName + " already in group " + groupName);
        exit(1);
    }


    // Create and write the attribute to the group

    H5::IntType integerType(H5::PredType::NATIVE_INT);
    H5::DataSpace attributeSpace(H5S_SCALAR);
    H5::Attribute attribute = group.createAttribute(attributeName.c_str(), integerType, attributeSpace);
    attribute.write(integerType, &attributeValue);
    attribute.close();

    // That's it

    attribute.close();
    group.close();

    return;
}













// HDF5File::writeAttribute()  for long integer attributes
//
// PURPOSE: add a long integer attribute to the HDF5 file.
//          The attribute is stored in "<GroupName>/<attributeName>".
//          
// INPUT: groupName:      string containing the full path of an existing group. Starts with "/".
//        attributeName:  string containing the name of the attribute
//        attributeValue: long integer containing the attribute value
//
// OUTPUT: None
//
// EXAMPLE: writeAttribute("/InputParameters/SeedParameters", "PhotonNoise", 1433237514)
//

void HDF5File::writeAttribute(string groupName, string attributeName, long attributeValue)
{
    // Complain if the file was not first opened
    
    if (!fileIsOpen)
    {
        Log.error("HDF5File::writeAttribute(): file is not open.");
        exit(1);
    }


    // Open the proper group where the input parameter belongs

    H5::Group group = file->openGroup(groupName.c_str());

    // Check whether the attribute is not already in the group. If so, complain.
    // The only way to do so, seems to try to access it and catch the exception
    // if it does not yet exist.

    bool attributeIsAlreadyInGroup = true;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the attribute

        H5::Attribute attr = group.openAttribute(attributeName.c_str());
    }
    catch (H5::AttributeIException error)
    {
        attributeIsAlreadyInGroup = false;
    }

    if (attributeIsAlreadyInGroup)
    {
        Log.error("HDF5File::writeAttribute(): attribute " + attributeName + " already in group " + groupName);
        exit(1);
    }


    // Create and write the attribute to the group

    H5::IntType integerType(H5::PredType::NATIVE_LONG);
    H5::DataSpace attributeSpace(H5S_SCALAR);
    H5::Attribute attribute = group.createAttribute(attributeName.c_str(), integerType, attributeSpace);
    attribute.write(integerType, &attributeValue);
    attribute.close();

    // That's it

    attribute.close();
    group.close();

    return;
}












// HDF5File::writeAttribute() for double-valued attributes
//
// PURPOSE: add a double-valued attribute to the HDF5 file.
//          The attribute is stored in "<GroupName>/<attributeName>".
//          
// INPUT: groupName:      string containing the full path of an existing group. Starts with "/".
//        attributeName:  string containing the name of the attribute
//        attributeValue: double containing the attribute value
//
// OUTPUT: None
//
// EXAMPLE: writeAttribute("/InputParameters/JitterParameters", "JitterYawRms", 0.01)
//

void HDF5File::writeAttribute(string groupName, string attributeName, double attributeValue)
{
    // Complain if the file was not first opened
    
    if (!fileIsOpen)
    {
        Log.error("HDF5File::writeAttribute(): file is not open.");
        exit(1);
    }

    // Open the proper group where the input parameter belongs

    H5::Group group = file->openGroup(groupName.c_str());

    // Check whether the attribute is not already in the group. If so, complain.
    // The only way to do so, seems to try to access it and catch the exception
    // if it does not yet exist.

    bool attributeIsAlreadyInGroup = true;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the attribute

        H5::Attribute attr = group.openAttribute(attributeName.c_str());
    }
    catch (H5::AttributeIException error)
    {
        attributeIsAlreadyInGroup = false;
    }

    if (attributeIsAlreadyInGroup)
    {
        Log.error("HDF5File::writeAttribute(): attribute " + attributeName + " already in group " + groupName);
        exit(1);
    }


    // Create and write the attribute to the group

    H5::IntType floatType(H5::PredType::NATIVE_DOUBLE);
    H5::DataSpace attributeSpace(H5S_SCALAR);
    H5::Attribute attribute = group.createAttribute(attributeName.c_str(), floatType, attributeSpace);
    attribute.write(floatType, &attributeValue);
    attribute.close();

    // That's it

    attribute.close();
    group.close();

    return;
}












// HDF5File::writeArray()  for 1D integer arrays
//
// PURPOSE: write a 1D array to a specified group in the HDF5 file.
//
// INPUT: groupName: full path of an existing HDF5 Group in the file. Starts with "/".
//        arrayName: unique name of the array in the group, e.g. "starIDs000001"
//        array:     1D integer native array
//        size:      number of elements in the array
//
// OUTPUT: None

void HDF5File::writeArray(string groupName, string arrayName, int* array, int size)
{
    // Create a DataSpace defining the shape and type of the data 

    unsigned int Ndimensions = 1;
    unsigned long long shape[Ndimensions];
    shape[0] = size;
    H5::DataSpace arraySpace(Ndimensions, shape);

    // Check if the array is not already in the file.
    // There seems to be only a dirty way of determining this:
    // try to access the dataset, and check if an exception is thrown.

    bool arrayIsAlreadyInFile = true;
    string arrayPath = groupName + "/" + arrayName;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the dataset

        H5::DataSet testDataset = file->openDataSet(arrayPath.c_str());
    }
    catch (H5::FileIException error)
    {
        arrayIsAlreadyInFile = false;
    }

    if (arrayIsAlreadyInFile)
    {
        Log.error("HDF5File::writeArray(): array " + groupName + "/" + arrayName + " already in file.");
        exit(1);
    }


    // Inside the Images group, make room for the image array

    H5::DataSet arrayDataset = file->createDataSet(arrayPath.c_str(), H5::PredType::NATIVE_INT, arraySpace);

    // Copy the data from our image into the HDF5 file

    arrayDataset.write(array, H5::PredType::NATIVE_INT);

    // That's it

    return;
}












// HDF5File::writeArray()  for 1D float arrays
//
// PURPOSE: write a 1D array to a specified group in the HDF5 file.
//
// INPUT: array: 1D float native array
//        groupName: name of an existing HDF5 Group in the file. Starts with "/".
//        arrayName: unique name of the array in the group, e.g. "starXcoordinates000001"
//
// OUTPUT: None

void HDF5File::writeArray(string groupName, string arrayName, float* array, int size)
{
     // Create a DataSpace defining the shape and type of the data 

    unsigned int Ndimensions = 1;
    unsigned long long shape[Ndimensions];
    shape[0] = size;
    H5::DataSpace arraySpace(Ndimensions, shape);

    // Check if the array is not already in the file.
    // There seems to be only a dirty way of determining this:
    // try to access the dataset, and check if an exception is thrown.

    bool arrayIsAlreadyInFile = true;
    string arrayPath = groupName + "/" + arrayName;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the dataset

        H5::DataSet testDataset = file->openDataSet(arrayPath.c_str());
    }
    catch (H5::FileIException error)
    {
        arrayIsAlreadyInFile = false;
    }

    if (arrayIsAlreadyInFile)
    {
        Log.error("HDF5File::writeArray(): array " + groupName + "/" + arrayName + " already in file.");
        exit(1);
    }


    // Inside the Images group, make room for the image array

    H5::DataSet arrayDataset = file->createDataSet(arrayPath.c_str(), H5::PredType::NATIVE_FLOAT, arraySpace);

    // Copy the data from our image into the HDF5 file

    arrayDataset.write(array, H5::PredType::NATIVE_FLOAT);

    // That's it

    return;
}











// HDF5File::writeArray()  for 1D double-valued arrays
//
// PURPOSE: write a 1D array to a specified group in the HDF5 file.
//
// INPUT: array: 1D double native array
//        groupName: name of an existing HDF5 Group in the file. Starts with "/".
//        arrayName: unique name of the array in the group, e.g. "DecOpticalAxis"
//
// OUTPUT: None

void HDF5File::writeArray(string groupName, string arrayName, double* array, int size)
{
    // Create a DataSpace defining the shape and type of the data 

    unsigned int Ndimensions = 1;
    unsigned long long shape[Ndimensions];
    shape[0] = size;
    H5::DataSpace arraySpace(Ndimensions, shape);

    // Check if the array is not already in the file.
    // There seems to be only a dirty way of determining this:
    // try to access the dataset, and check if an exception is thrown.

    bool arrayIsAlreadyInFile = true;
    string arrayPath = groupName + "/" + arrayName;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the dataset

        H5::DataSet testDataset = file->openDataSet(arrayPath.c_str());
    }
    catch (H5::FileIException error)
    {
        arrayIsAlreadyInFile = false;
    }

    if (arrayIsAlreadyInFile)
    {
        Log.error("HDF5File::writeArray(): array " + groupName + "/" + arrayName + " already in file.");
        exit(1);
    }


    // Inside the Images group, make room for the image array

    H5::DataSet arrayDataset = file->createDataSet(arrayPath.c_str(), H5::PredType::NATIVE_DOUBLE, arraySpace);

    // Copy the data from our image into the HDF5 file

    arrayDataset.write(array, H5::PredType::NATIVE_DOUBLE);

    // That's it

    return;
}











// HDF5File::writeArray()  for 2D float armadillo arrays
//
// PURPOSE: write a 2D armadillo array to a specified group in the HDF5 file.
//
// INPUT: array: 2D armadillo array
//        groupName: name of an existing HDF5 Group in the file. Starts with "/".
//        arrayName: unique name of the array in the group, e.g. "image000001"
//
// OUTPUT: None

void HDF5File::writeArray(string groupName, string arrayName, arma::Mat<float>& A)
{
    // Sanity check on the shape of the array

    if ((A.n_rows == 0) && (A.n_cols == 0))
    {
        Log.error("HDF5File::writeArray(): encountered array with shape (0,0)");
        exit(1);
    }

    // Create a DataSpace defining the shape and type of the data 

    unsigned int Ndimensions = 2;
    unsigned long long shape[Ndimensions];
    shape[0] = A.n_rows;
    shape[1] = A.n_cols;
    H5::DataSpace arraySpace(Ndimensions, shape);

    // Check if the array is not already in the file.
    // There seems to be only a dirty way of determining this:
    // try to access the dataset, and check if an exception is thrown.

    bool arrayIsAlreadyInFile = true;
    string arrayPath = groupName + "/" + arrayName;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the dataset

        H5::DataSet testDataset = file->openDataSet(arrayPath.c_str());
    }
    catch (H5::FileIException error)
    {
        arrayIsAlreadyInFile = false;
    }

    if (arrayIsAlreadyInFile)
    {
        Log.error("HDF5File::writeArray(): array " + groupName + "/" + arrayName + " already in file.");
        exit(1);
    }


    // Inside the Images group, make room for the image array

    H5::DataSet arrayDataset = file->createDataSet(arrayPath.c_str(), H5::PredType::NATIVE_FLOAT, arraySpace);

    // Copy the data from our image into the HDF5 file

    arrayDataset.write(A.memptr(), H5::PredType::NATIVE_FLOAT);


    // That's it

    return;
}










// HDF5File::readArray() for 2D float armadillo arrays
//
// PURPOSE: read a 2D array from a specified group in the HDF5 file into an
//          armadillo array.
//
// INPUT: array: 2D armadillo array. Previous contents will be lost.
//        groupName: name of an existing HDF5 Group in the file. Starts with "/".
//        arrayName: unique name of the array in the group, e.g. "image000001"
//
// OUTPUT: None

void HDF5File::readArray(string groupName, string arrayName, arma::Mat<float>& A)
{
    // Construct the path of the dataset in the HDF5 file

    string arrayPath = groupName + "/" + arrayName;

    // Try to open the dataset

    H5::DataSet dataset;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the dataset

        dataset = file->openDataSet(arrayPath.c_str());
    }
    catch (H5::FileIException error)
    {
        Log.error("readArray(): " + arrayPath + " not in file.");
        exit(1);
    }

    // Find out the number of rows and columns of the dataset

    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t shape[2];
    unsigned int Ndimensions = dataspace.getSimpleExtentDims(shape, NULL);
    int Nrows = shape[0];
    int Ncolumns = shape[1];

    // Reset the size for the 2D array

    A.reset();
    A.set_size(Nrows, Ncolumns);

    // Read the HDF5 dataset into the array

    dataset.read(A.memptr(), H5::PredType::NATIVE_FLOAT);

    // That's it

    return;
}

