/**
 * @class      ConfigurationParameters 
 * 
 * @brief      Parse the input file and make all input parameters available to the Simulator.
 * 
 * @details
 * 
 * The ConfigurationParameters provide an easy way to load and use input parameters
 * in the simulator. All parameters are loaded from a single YAML input file. The file is loaded 
 * and parsed by the constructor. 
 */
#include <string>
#include <list>

#include "FileUtilities.h"
#include "StringUtilities.h"
#include "Exceptions.h"
#include "ConfigurationParameters.h"
#include "Logger.h"

using namespace std;


ConfigurationParameters::ConfigurationParameters() {}





/**
 * @brief      Loads the input file. The input is expected to be a YAML file.
 *
 * @exception  IOException is thrown when the file does not exist
 * 
 * @param[in]  name  Filename of the input file for PlatoSim3
 */
ConfigurationParameters::ConfigurationParameters(const char* name) : ConfigurationParameters::ConfigurationParameters(string(name)) {}





/**
 * @brief      Loads the input file. The input is expected to be a YAML file.
 *
 * @exception  IOException is thrown when the file does not exist
 * 
 * @param[in]  name  Filename of the input file for PlatoSim3
 */
ConfigurationParameters::ConfigurationParameters(const string &name)
{
    if ( ! FileUtilities::fileExists(name) )
    {
        Log.warning("Error: Filename \"" + name + "\" does not exist.");
        throw IOException("File passed as an argument to ConfigurationParameters does not exist.");
    }

    filename = name;
    config = YAML::LoadFile(name);

}





/**
 * @brief      Return the boolean value for the specified parameter.
 *
 * @details    
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 * is part of the group ObservingParameters, then the key to get the value for this
 * parameter would be "ObservingParameters/ExposureTime".
 * 
 * @param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * @returns A boolean value for the given parameter
 */
bool ConfigurationParameters::getBoolean(const string &key)
{
    vector<string> fields = StringUtilities::split(key, '/');

    if (fields.size() > 1)
    {
        YAML::Node node = config[fields[0]];
        return node[fields[1]].as<bool>();
    }
    else 
    {
        return config[key].as<bool>();
    }

}





/**
 * @brief      Return the integer value for the specified parameter.
 *
 * @details    
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 * is part of the group ObservingParameters, then the key to get the value for this
 * parameter would be "ObservingParameters/ExposureTime".
 * 
 * @param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * @returns    An integer value for the given parameter
 */
int ConfigurationParameters::getInteger(const string &key)
{
    vector<string> fields = StringUtilities::split(key, '/');

    if (fields.size() > 1)
    {
        YAML::Node node = config[fields[0]];
        return node[fields[1]].as<int>();
    }
    else 
    {
        return config[key].as<int>();
    }

}





/**
 * @brief      Return the double value for the specified parameter.
 *
 * @details    
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 * is part of the group ObservingParameters, then the key to get the value for this
 * parameter would be "ObservingParameters/ExposureTime".
 * 
 * @param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * @returns    A double value for the given parameter
 */
double ConfigurationParameters::getDouble(const string &key)
{
    vector<string> fields = StringUtilities::StringUtilities::StringUtilities::split(key, '/');

    if (fields.size() > 1)
    {
        YAML::Node node = config[fields[0]];
        return node[fields[1]].as<double>();
    }
    else 
    {
        return config[key].as<double>();
    }

}





/**
 * @brief      Return the string value for the specified parameter.
 *
 * @details    
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 * is part of the group ObservingParameters, then the key to get the value for this
 * parameter would be "ObservingParameters/ExposureTime".
 * 
 * @param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * @returns    A string value for the given parameter
 */
string ConfigurationParameters::getString(const string &key) 
{
    vector<string> fields = StringUtilities::StringUtilities::split(key, '/');

    if (fields.size() > 1)
    {
        YAML::Node node = config[fields[0]];
        return node[fields[1]].as<string>();
    }
    else 
    {
        return config[key].as<string>();
    }

}





/**
 * @brief      Return the absolute filename for the given parameter.
 * 
 * @details    
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 * is part of the group ObservingParameters, then the key to get the value for this
 * parameter would be "ObservingParameters/ExposureTime".
 * 
 * When the parameter contains an absolute filename, that value is returned. 
 * If the parameter contains a relative path, the filename is preceeded by the value of the 
 * General/ProjectLocation parameter. An absolute path starts with a '/' character, 
 * otherwise the path is considered relative.
 *
 * @param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * @returns    An absolute filename
 */
string ConfigurationParameters::getAbsoluteFileName(const string &key) 
{
    string filename;
    vector<string> fields = StringUtilities::split(key, '/');

    if (fields.size() > 1)
    {
        YAML::Node node = config[fields[0]];
        filename = node[fields[1]].as<string>();
    }
    else 
    {
        filename = config[key].as<string>();
    }

    if (FileUtilities::isRelative(filename))
    {
        string projectLocation = this->getString("General/ProjectLocation");
        return projectLocation + "/" + filename;
    }

    return filename;
}




ConfigurationParameters::~ConfigurationParameters() {}





