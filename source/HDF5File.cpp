/**
 * \class HDF5File
 * 
 * \brief      Provides an application level interface to the HDF5 C++ API
 * 
 * \details
 * 
 * This class provides a convenient application level interface to the HDF5 C++ wrapper to 
 * the HDF C library that is developed by the HDFGroup.
 * 
 * Use this class to access (read/write) HDF5 files from your code. Do not use the HDF5 C++
 * wrapper.
 */
 
#include "HDF5File.h"



/**
 * \brief      Default Constructor
 */
 
HDF5File::HDF5File()
: file(NULL), fileIsOpen(false)
{

}







/**
 * \brief  Open file constructor
 *
 * \param  filename   Absolute path of the file
 * \param  readonly   True if an existing file should only be read and not written.
 *                    False otherwise. Ignored if the file does not exist yet.
 */

HDF5File::HDF5File(string filename, bool readonly)
: file(NULL), fileIsOpen(false)
{
    H5::Exception::dontPrint();
    open(filename, readonly);
}







/**
 * \brief      Destructor
 */

HDF5File::~HDF5File()
{
    if (fileIsOpen)
    {
        close();
    }
    delete file;
}








/**
 * \brief Open an HDF5 File. If it doesn't already exist, create it.
 * 
 * 
 * \param filename   Absolute path of the file.
 * \param readonly   True if an existing file should only be read and not written.
 *                   False otherwise. Ignored if the file does not exist yet.
 */

void HDF5File::open(string filename, bool readonly)
{
    // If there was already a file opened, complain and exit

    if (fileIsOpen)
    {
        close();
        throw H5FileException("HDF5File::open(): file " + filename + " is already open.");
    }

    // Check if the file already exists.
    // If so, use append mode, if not

    if (fileExists(filename))
    {
        if (readonly)
        {
            // Open existing HDF5 file to read.

            file = new H5::H5File(filename.c_str(), H5F_ACC_RDONLY);

            Log.info("HDF5File: opened existing HDF5 file " + filename + " to read");

        }
        else
        {
            // Open an existing HDF5 file to both read and write
    
            file = new H5::H5File(filename.c_str(), H5F_ACC_RDWR);

            Log.info("HDF5File: opened existing HDF5 file " + filename + " to read/write");            
        }
    }
    else
    {
        // The file doesn't exist yet. Create a new one.
        // New files are always opened for read/write.

        file = new H5::H5File(filename.c_str(), H5F_ACC_TRUNC);

        Log.info("HDF5File: opened new HDF5 file " + filename + " to read/write");
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
        file = nullptr;
        fileIsOpen = false;
    }
}











/**
 * \brief      Check if the HDF5 file has a group with the given name
 *
 * \param[in]  groupName    the full name of the group
 *
 * \return     true if the group exists in this file, false otherwise
 */

bool HDF5File::hasGroup(string groupName)
{
    try
    {
        H5::Group parentGroup = file->openGroup(groupName.c_str());
        return true;
    }
    catch (H5::Exception ex) {}

    return false;
}






/**
 * \brief      Check if the given group has a dataset with the given name
 *
 * \param[in]  groupName    the full name of the group that contains the dataset
 * \param[in]  datasetName  the name of the dataset
 *
 * \return     true if the dataset exists in this group, false otherwise
 * 
 * \exception  H5GroupException thrown when the group is unknown to the HDF5 file
 */
 
bool HDF5File::hasDataset(string groupName, string datasetName)
{
    if (hasGroup(groupName))
    {
        H5::Group group = file->openGroup(groupName.c_str());

        try
        {
            H5::DataSet dataset = group.openDataSet(datasetName.c_str());
            return true;
        }
        catch (H5::Exception ex) {}

        return false;
    }
    else
    {
        throw H5GroupException("HDF5File: Unknown group (" + groupName + ") in HDF5 file " + file->getFileName());
    }

    return false;
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
        throw H5FileException("HDF5File::writeAttribute(): file " + file->getFileName() + " is not open.");
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
        throw H5AttributeException("HDF5File::writeAttribute(): attribute " + attributeName + " already in group " + groupName);
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
        throw H5FileException("HDF5File::writeAttribute(): file " + file->getFileName() + " is not open.");
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
        throw H5AttributeException("HDF5File::writeAttribute(): attribute " + attributeName + " already in group " + groupName);
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
        throw H5FileException("HDF5File::writeAttribute(): file " + file-> getFileName() + " is not open.");
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
        throw H5AttributeException("HDF5File::writeAttribute(): attribute " + attributeName + " already in group " + groupName);
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
        throw H5FileException("HDF5File::writeAttribute(): file " + file->getFileName() + " is not open.");
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
        throw H5AttributeException("HDF5File::writeAttribute(): attribute " + attributeName + " already in group " + groupName);
    }


    // Create and write the attribute to the group

    H5::FloatType floatType(H5::PredType::NATIVE_DOUBLE);
    H5::DataSpace attributeSpace(H5S_SCALAR);
    H5::Attribute attribute = group.createAttribute(attributeName.c_str(), floatType, attributeSpace);
    attribute.write(floatType, &attributeValue);
    attribute.close();

    // That's it

    attribute.close();
    group.close();

    return;
}










/**
 * \brief Add a boolean attribute to the HDF5 file. 
 *        The attribute is stored in "<GroupName>/<attributeName>".
 *        
 * \param groupName        String containing the full path of an existing group. Starts with "/".
 * \param attributeName    String containing the name of the attribute
 * \param attributeValue   Integer containing the attribute value
 * 
 * \example  writeAttribute("/InputParameters/Platform", "UseJitter", true) 
 * 
 * \exception  H5FileException   If the HDF5 file has not been opened
 * \exception  H5FileException   If the attribute already existed in the HDF5 file
 */

void HDF5File::writeAttribute(string groupName, string attributeName, bool attributeValue)
{
    // Complain if the file was not first opened
    
    if (!fileIsOpen)
    {
        string errorMessage = "HDF5File::writeAttribute(): file is not open.";
        Log.error(errorMessage);
        throw H5FileException(errorMessage);
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
        string errorMessage = "HDF5File::writeAttribute(): attribute " + attributeName + " already in group " + groupName;
        Log.error(errorMessage);
        throw H5FileException(errorMessage);
    }


    // Create and write the attribute to the group
    // First convert from boolean (true/false) to integer (1/0)

    int value;
    if (attributeValue)
    {
        value = 1;
    }
    else
    {
        value = 0;
    }

    H5::IntType integerType(H5::PredType::NATIVE_INT);
    H5::DataSpace attributeSpace(H5S_SCALAR);
    H5::Attribute attribute = group.createAttribute(attributeName.c_str(), integerType, attributeSpace);
    attribute.write(integerType, &value);
    attribute.close();

    // That's it

    attribute.close();
    group.close();

    return;
}













// HDF5File::writeAttribute() for an attribute of vector <double>
//
// PURPOSE: add an attribute of type vector <double> to the HDF5 file.
//          The attribute is stored in "<GroupName>/<attributeName>".
//          
// INPUT: groupName:      string containing the full path of an existing group. Starts with "/".
//        attributeName:  string containing the name of the attribute
//        attributeValue: vector <double> containing the attribute values
//
// OUTPUT: None
//
// EXAMPLE: writeAttribute("/InputParameters/Camera/FieldDistortion", "Coefficients", values)
//

void HDF5File::writeAttribute(string groupName, string attributeName, vector<double> attributeValue)
{
    // Complain if the file was not first opened
    
    if (!fileIsOpen)
    {
        throw H5FileException("HDF5File::writeAttribute(): file " + file->getFileName() + " is not open.");
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
        throw H5AttributeException("HDF5File::writeAttribute(): attribute " + attributeName + " already in group " + groupName);
    }


    // Create and write the attribute to the group. The attribute is a vector <double> so we need
    // to specify the rank and the dimension of the vector. Since a vector is guaranteed to store 
    // their elements contiguously, we can just pass the pointer to the first element when writing 
    // the attribute.

    try
    {
        const int rank = 1;
        hsize_t dimension[1];
        dimension[0] = attributeValue.size();
        H5::FloatType floatType(H5::PredType::NATIVE_DOUBLE);
        H5::DataSpace attributeSpace(rank, dimension);
        H5::Attribute attribute = group.createAttribute(attributeName.c_str(), floatType, attributeSpace);
        attribute.write(floatType, &attributeValue[0]);
        attribute.close();
    }
    catch(H5::AttributeIException error)
    {
        throw H5AttributeException("HDF5File::writeAttribute(): attribute " + attributeName + " write error " + error.getCDetailMsg());
    }

    // That's it

    group.close();

    return;
}














/**
 * \brief      read a double-valued attribute that is associated with groupName
 *
 * \param[in]  groupName      string containing the full path of an existing group. Starts with "/".
 * \param[in]  attributeName  string containing the name of the attribute
 *
 * \return     the value of the attribute
 * 
 * \exception  H5FileException      if the HDF5 file has not been opened
 * \exception  H5GroupException     if the group is unknown to the HDF5 file
 * \exception  H5AttributeException if there is no attribute with the given name
 */

double HDF5File::readDoubleGroupAttribute(string groupName, string attributeName)
{
    // Complain if the file was not first opened
    
    if (!fileIsOpen)
    {
        throw H5FileException("HDF5File: The file (" + file->getFileName() + ") has not been opened.");
    }

    // Open the proper group where the attribute is associated

    H5::Group group;
    if (hasGroup(groupName))
    {
        group = file->openGroup(groupName.c_str());
    }
    else 
    {
        throw H5GroupException("HDF5File: Unknown group (" + groupName + ") in HDF5 file " + file->getFileName());
    }

    // Check whether the attribute is in the group by trying to read it.
    // If not, raise an exception.

    H5::Attribute attr;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the attribute

        attr = group.openAttribute(attributeName.c_str());
    }
    catch (H5::AttributeIException error)
    {
        throw H5AttributeException("HDF5File: Unknown Attribute (" + attributeName + ") in the group " + groupName + " for HDF5 file " + file->getFileName());
    }

    double value = 0.0;

    H5::DataType type = attr.getDataType();
    attr.read(type, &value);    

    return value;
}












/**
 * \brief      read a integer-valued attribute that is associated with groupName
 *
 * \param[in]  groupName      string containing the full path of an existing group. Starts with "/".
 * \param[in]  attributeName  string containing the name of the attribute
 *
 * \return     the value of the attribute
 * 
 * \exception  H5FileException      if the HDF5 file has not been opened
 * \exception  H5GroupException     if the group is unknown to the HDF5 file
 * \exception  H5AttributeException if there is no attribute with the given name
 */

int HDF5File::readIntegerGroupAttribute(string groupName, string attributeName)
{
    // Complain if the file was not first opened
    
    if (!fileIsOpen)
    {
        throw H5FileException("HDF5File: The file (" + file->getFileName() + ") has not been opened.");
    }

    // Open the proper group where the attribute is associated

    H5::Group group;
    if (hasGroup(groupName))
    {
        group = file->openGroup(groupName.c_str());
    }
    else 
    {
        throw H5GroupException("HDF5File: Unknown group (" + groupName + ") in HDF5 file " + file->getFileName());
    }

    // Check whether the attribute is in the group by trying to read it.
    // If not, raise an exception.

    H5::Attribute attr;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the attribute

        attr = group.openAttribute(attributeName.c_str());
    }
    catch (H5::AttributeIException error)
    {
        throw H5AttributeException("HDF5File: Unknown Attribute (" + attributeName + ") in the group " + groupName + " for HDF5 file " + file->getFileName());
    }

    int value = 0.0;

    H5::DataType type = attr.getDataType();
    attr.read(type, &value);    

    return value;
}














/**
 * \brief      read a double-valued attribute that is associated with dataset
 *
 * \param[in]  groupName      string containing the full path of an existing group. Starts with "/".
 * \param[in]  datasetName    string containing the name of the dataSet where the attribute is attached
 * \param[in]  attributeName  string containing the name of the attribute
 *
 * \return     the value of the attribute
 * 
 * \exception  H5FileException      if the HDF5 file has not been opened
 * \exception  H5GroupException     if the group is unknown to the HDF5 file
 * \exception  H5DatasetException   if the dataset is not known to the group
 * \exception  H5AttributeException if there is no attribute with the given name
 */

double HDF5File::readDoubleDatasetAttribute(string groupName, string datasetName, string attributeName)
{
    // Complain if the file was not first opened
    
    if (!fileIsOpen)
    {
        throw H5FileException("HDF5File: The file (" + file->getFileName() + ") has not been opened.");
    }

    // Open the proper group where the attribute is associated

    H5::Group group;
    if (hasGroup(groupName))
    {
        group = file->openGroup(groupName.c_str());
    }
    else 
    {
        throw H5GroupException("HDF5File: Unknown group (" + groupName + ") in HDF5 file " + file->getFileName());
    }

    H5::DataSet dataset;
    if (hasDataset(groupName, datasetName))
    {
        dataset = group.openDataSet(datasetName);
    }
    else 
    {
        throw H5DatasetException("HDF5File: Unknown dataset (" + datasetName + ") in group (" + groupName + ") in HDF5 file " + file->getFileName());
    }

    // Check whether the attribute is in the group by trying to read it.
    // If not, raise an exception.

    H5::Attribute attr;

    try 
    {  
        // Turn off the auto-printing when an exception is raised

        H5::Exception::dontPrint();

        // Try to open the attribute

        attr = dataset.openAttribute(attributeName.c_str());
    }
    catch (H5::AttributeIException error)
    {
        throw H5AttributeException("HDF5File: Unknown Attribute (" + attributeName + ") in the group " + groupName + " for HDF5 file " + file->getFileName());
    }

    double value = 0.0;

    H5::DataType type = attr.getDataType();
    attr.read(type, &value);

    return value;
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
        throw H5GroupException("HDF5File::writeArray(): array " + groupName + "/" + arrayName + " already in file.");
    }


    // Inside the Images group, make room for the image array

    H5::DataSet arrayDataset = file->createDataSet(arrayPath.c_str(), H5::PredType::NATIVE_INT, arraySpace);

    // Copy the data from our image into the HDF5 file

    arrayDataset.write(array, H5::PredType::NATIVE_INT);

    // That's it

    return;
}













/**
 * \brief Write a 1D unsigned integer array to a specified group in the HDF5 file
 * 
 * \param groupName  Full path of an existing HDF5 Group in the file. Starts with "/". 
 * \param arrayName  Unique name of the array in the group 
 * \param array      1d unsigned integer native array
 * \param size       Number of elements in the array
 */

void HDF5File::writeArray(string groupName, string arrayName, unsigned int* array, int size)
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
        string errorMessage = "HDF5File::writeArray(): array " + groupName + "/" + arrayName + " already in file.";
        Log.error(errorMessage);
        throw H5FileException(errorMessage);
    }


    // Inside the Images group, make room for the image array

    H5::DataSet arrayDataset = file->createDataSet(arrayPath.c_str(), H5::PredType::NATIVE_UINT, arraySpace);

    // Copy the data from our image into the HDF5 file

    arrayDataset.write(array, H5::PredType::NATIVE_UINT);

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
        throw H5GroupException("HDF5File::writeArray(): array " + groupName + "/" + arrayName + " already in file.");
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
        throw H5GroupException("HDF5File::writeArray(): array " + groupName + "/" + arrayName + " already in file.");
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
        throw H5FileException("HDF5File::writeArray(): encountered array with shape (0,0)");
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
        throw H5GroupException("HDF5File::writeArray(): array " + groupName + "/" + arrayName + " already in file.");
    }

    // Inside the Images group, make room for the image array

    H5::DataSet arrayDataset = file->createDataSet(arrayPath.c_str(), H5::PredType::NATIVE_FLOAT, arraySpace);
    
    // Copy the Armadillo array to a vector, because the internally Armadillo stores the data column-major
    // while HDF5 assumes data to be stored row-major

    vector<float> temp(A.n_rows * A.n_cols);
    for (int n = 0; n < A.n_rows; n++)
    {
        for(int k = 0; k < A.n_cols; k++)
        {
            const int nk = n * A.n_cols + k;
            temp[nk] = A(n,k);
        }
    }

    // Copy the data from our image into the HDF5 file

    arrayDataset.write(temp.data(), H5::PredType::NATIVE_FLOAT);

    // That's it

    return;
}










// HDF5File::readArray() for 2D float armadillo arrays
//
// PURPOSE: read a 2D array from a specified group in the HDF5 file into an
//          armadillo array.
//
// INPUT: array:     2D armadillo array. Previous contents will be lost.
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
        throw H5DatasetException("HDF5File::readArray(): " + arrayPath + " not in file.");
    }

    // Find out the number of rows and columns of the dataset

    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t shape[2];
    unsigned int Ndimensions = dataspace.getSimpleExtentDims(shape, NULL);
    int Nrows = shape[0];
    int Ncolumns = shape[1];


    // Make a temporary vector<> object. We cannot read directly into the Armadillo array
    // because the latter are column-major while HDF5 stores data row-major

    vector<float> temp(Nrows*Ncolumns);

    // Read the HDF5 dataset into the array

    dataset.read(temp.data(), H5::PredType::NATIVE_FLOAT);

    // Reset the size for the 2D Armadillo array, and copy the data.

    A.reset();
    A.set_size(Nrows, Ncolumns);

    for (int n = 0; n < A.n_rows; n++)
    {
        for(int k = 0; k < A.n_cols; k++)
        {
            const int nk = n * A.n_cols + k;
            A(n,k) = temp[nk];                   // Is stored column-major
        }
    }

    // That's it

    return;
}










/**
 * \brief  Read a 1D double array from a specified group in the HDF5 file into a vector<double>.
 * 
 * \param groupName  Name of an existing HDF5 Group in the file. Starts with "/".
 * \param arrayName  Unique name of the array in the group, e.g. "skyBackground"
 * \param vec        C++ vector<double>. Previous contents will be lost.
 * 
 */

void HDF5File::readArray(string groupName, string arrayName, vector<double> &vec)
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
        throw H5DatasetException("HDF5File::readArray(): " + arrayPath + " not in file.");
    }

    // Find out the size of the dataset

    H5::DataSpace dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1)
    {
        throw H5DatasetException("HDF5File::readArray(): " + arrayPath + " is not 1D.");
    }

    hsize_t shape[1];
    unsigned int Ndimensions = dataspace.getSimpleExtentDims(shape, NULL);
    int size = shape[0];


    // Ensure that the vector has enough space to read all data

    vec.clear();
    vec.resize(size);

    // Read the HDF5 dataset into the array

    dataset.read(vec.data(), H5::PredType::NATIVE_DOUBLE);

    // That's it

    return;
}












/**
 * \brief  Read a 1D (unsigned int) array from a specified group in the HDF5 file into a vector<unsigned int>
 * 
 * \param groupName  Name of an existing HDF5 Group in the file. Starts with "/".
 * \param arrayName  Unique name of the array in the group, e.g. "skyBackground"
 * \param vec        C++ vector<double>. Previous contents will be lost.
 * 
 */

void HDF5File::readArray(string groupName, string arrayName, vector<unsigned int> &vec)
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
        throw H5DatasetException("HDF5File::readArray(): " + arrayPath + " not in file.");
    }

    // Find out the size of the dataset

    H5::DataSpace dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();
    if (rank != 1)
    {
        throw H5DatasetException("HDF5File::readArray(): " + arrayPath + " is not 1D.");
    }

    hsize_t shape[1];
    unsigned int Ndimensions = dataspace.getSimpleExtentDims(shape, NULL);
    int size = shape[0];


    // Ensure that the vector has enough space to read all data

    vec.clear();
    vec.resize(size);

    // Read the HDF5 dataset into the array

    dataset.read(vec.data(), H5::PredType::NATIVE_UINT);

    // That's it

    return;
}














// fileExists()
//
// PURPOSE: Check if file exists. This is done by trying to open the file,
//          and verifying if the "good" flag is raised.
//
// INPUT: filename: full path of the file
//
// OUTPUT: true if the file exists, false otherwise

bool fileExists(string filename)
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
