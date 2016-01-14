#include <string>
#include <list>

#include "FileUtilities.h"
#include "StringUtilities.h"
#include "Exceptions.h"
#include "JSonFormat.h"
#include "YamlFormat.h"
#include "ConfigurationParameters.h"
#include "Logger.h"



ConfigurationParameters::ConfigurationParameters() {}

ConfigurationParameters::ConfigurationParameters(const char* name) : ConfigurationParameters::ConfigurationParameters(std::string(name)) {}

ConfigurationParameters::ConfigurationParameters(const std::string &name)
{
    if ( ! fileExists(name))
    {
        Log.warning("Error: Filename \"" + name + "\" does not exist.");
        throw IOException("File passed as an argument to ConfigurationParameters does not exist.");
    }

    filename = name;

    if (ends_with(filename, ".json"))
    {
        format = new JSonFormat(filename);
    }
    else if (ends_with(filename, ".yaml"))
    {
        format = new YamlFormat(filename);
    }

}





int ConfigurationParameters::getInteger(const std::string &key)
{
    return format->getInteger(key);
}





double ConfigurationParameters::getDouble(const std::string &key)
{
    return format->getDouble(key);
}





/**
 * PURPOSE: Return the string value for the specified parameter
 *
 * INPUTS:  parameterName The name of a parameter used in the PLATO Simulator
 *
 * OUTPUTS: A string value for the given parameter
 */
std::string ConfigurationParameters::getString(const std::string &key) 
{
    return format->getString(key);
}





ConfigurationParameters::~ConfigurationParameters() {}





