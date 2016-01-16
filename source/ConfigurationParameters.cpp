#include <string>
#include <list>

#include "FileUtilities.h"
#include "StringUtilities.h"
#include "Exceptions.h"
#include "ConfigurationParameters.h"
#include "Logger.h"

using namespace std;



ConfigurationParameters::ConfigurationParameters() {}

ConfigurationParameters::ConfigurationParameters(const char* name) : ConfigurationParameters::ConfigurationParameters(string(name)) {}

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
 * PURPOSE: Return the boolean value for the specified parameter.
 *
 * INPUTS:  key: The name of a parameter used in the PLATO Simulator
 *
 * OUTPUTS: A boolean value for the given parameter
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
 * PURPOSE: Return the integer value for the specified parameter.
 *
 * INPUTS:  key: The name of a parameter used in the PLATO Simulator
 *
 * OUTPUTS: An integer value for the given parameter
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
 * PURPOSE: Return the double value for the specified parameter.
 *
 * INPUTS:  key: The name of a parameter used in the PLATO Simulator
 *
 * OUTPUTS: A double value for the given parameter
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
 * PURPOSE: Return the string value for the specified parameter.
 *
 * INPUTS:  key: The name of a parameter used in the PLATO Simulator
 *
 * OUTPUTS: A string value for the given parameter
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
 * PURPOSE: Return the absolute filename for the given parameter. When the parameter contains an
 *          absolute filename, that value is returned. If the parameter contains a relative path,
 *          the filename is constructed from the General/ProjectLocation and the relative filename
 *          of the parameter.
 *
 * INPUTS:  key: The name of a parameter used in the PLATO Simulator
 *
 * OUTPUTS: An absolute filename
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





