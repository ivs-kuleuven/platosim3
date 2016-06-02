/**
 * \class      ConfigurationParameters 
 * 
 * \brief      Parse the input file and make all input parameters available to the Simulator.
 * 
 * \details
 * 
 * The ConfigurationParameters provide an easy way to load and use input parameters
 * in the simulator. All parameters are loaded from a single YAML input file. The file is loaded 
 * and parsed by the constructor. 
 * 
 * TODO: Add explanantion of the keys, i.e. that they need to contain the group or section e.g. "General/ProjectLocation".
 * 
 */
#include "ConfigurationParameters.h"







/**
 * \brief      Default constructor
 */
ConfigurationParameters::ConfigurationParameters() {}









/**
 * \brief      Constructor. Loads the input file. The input is expected to be a YAML file.
 *
 * \exception  IllegalArgumentException is thrown when the file does not exist
 * 
 * \param[in]  name  Filename of the input file for PlatoSim3
 */
ConfigurationParameters::ConfigurationParameters(const string &name)
{
    if ( ! FileUtilities::fileExists(name) )
    {
        throw IllegalArgumentException("ConfigurationParameters: File (" + name + ") passed as an argument to ConfigurationParameters does not exist.");
    }

    filename = name;
    config = YAML::LoadFile(name);

}





/**
 * \brief      Destructor
 */
ConfigurationParameters::~ConfigurationParameters()
{
    // TODO: Find out if it is needed to close the YAML file and how to do this.
}






/**
 * \brief      Return the boolean value for the specified parameter.
 *
 * \details    
 * 
 * A boolean value is always parsed from a string that can have the following values: "yes"/"no", "y"/"n", "true"/"false", "on"/"off".
 * All the previous values can be Capitalized (i.e. start with an upper case letter) or can be all caps.
 * 
 * A boolean value can also be parsed from the integers 1 or 0.
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 * is part of the group ObservingParameters, then the key to get the value for this
 * parameter would be "ObservingParameters/ExposureTime".
 * 
 * \param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * \returns A boolean value for the given parameter
 */
bool ConfigurationParameters::getBoolean(const string &key)
{
    bool value;

    YAML::Node node = getNode(key);

    // First try to convert to a boolean, values should be true/false, on/off, yes/no, ..
    // Trying to convert from an integer into a boolean will fail with an Exception.

    try
    {
        value = node.as<bool>();
    }
    catch(YAML::Exception ex)
    {
        // Try to convert from an integer 0 or 1 

        try
        {
            int iValue = node.as<int>();
            if (iValue == 0)
                value = false;
            else if (iValue == 1)
                value = true;
            else
                throw ConfigurationException("ConfigurationParameters: expected boolean value for key=\"" + key + "\", while value=" + node.as<string>());
        }
        catch(YAML::Exception ex)
        {
            throw ConfigurationException("ConfigurationParameters: cannot convert key to boolean: key=\"" + key + "\", value=" + node.as<string>());
        }

    }

    return value;
}





/**
 * \brief      Return the integer value for the specified parameter.
 *
 * \details    
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 * is part of the group ObservingParameters, then the key to get the value for this
 * parameter would be "ObservingParameters/ExposureTime".
 * 
 * \param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * \returns    An integer value for the given parameter
 */
int ConfigurationParameters::getInteger(const string &key)
{
    int value;

    YAML::Node node = getNode(key);

    try
    {
        value = node.as<int>();
    }
    catch(YAML::Exception ex)
    {
        throw ConfigurationException("ConfigurationParameters: cannot convert key to integer: key=\"" + key + "\", value= " + node.as<string>());
    }

    return value;
}












/**
 * \brief      Return the long integer value for the specified parameter.
 *
 * \details    
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter
 * fullWellSaturationLimit is part of the group CCD, then the key to get the value 
 * for this parameter would be "CCD/fullWellSaturationLimit".
 * 
 * \param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * \returns    An integer value for the given parameter
 */
long ConfigurationParameters::getLong(const string &key)
{
    long value;

    YAML::Node node = getNode(key);

    try
    {
        value = node.as<long>();
    }
    catch(YAML::Exception ex)
    {
        throw ConfigurationException("ConfigurationParameters: cannot convert key to integer type long: key=\"" + key + "\", value= " + node.as<string>());
    }

    return value;
}














/**
 * \brief      Return the double value for the specified parameter.
 *
 * \details    
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 * is part of the group ObservingParameters, then the key to get the value for this
 * parameter would be "ObservingParameters/ExposureTime".
 * 
 * \param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * \returns    A double value for the given parameter
 */
double ConfigurationParameters::getDouble(const string &key)
{
    double value;

    YAML::Node node = getNode(key);

    try
    {
        value = node.as<double>();
    }
    catch(YAML::Exception ex)
    {
        throw ConfigurationException("ConfigurationParameters: cannot convert key to double: key=\"" + key + "\", value= " + node.as<string>());
    }

    return value;
}





/**
 * \brief      Return a vector of doubles for the specified parameter.
 *
 * \details    
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 * is part of the group ObservingParameters, then the key to get the value for this
 * parameter would be "ObservingParameters/ExposureTime".
 * 
 * \param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * \returns    A vector of double values for the given parameter
 */
vector <double> ConfigurationParameters::getDoubleVector(const string &key)
{

    YAML::Node valuesNode = getNode(key);

    unsigned short nValues = valuesNode.size();

    vector<double> values(nValues);

    for (unsigned short idx = 0; idx < nValues; ++idx){
            values[idx] = valuesNode[idx].as<double>();
    }

    return values;
}




/**
 * \brief      Return the string value for the specified parameter.
 *
 * \details    
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 * is part of the group ObservingParameters, then the key to get the value for this
 * parameter would be "ObservingParameters/ExposureTime".
 * 
 * \param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * \returns    A string value for the given parameter
 */
string ConfigurationParameters::getString(const string &key) 
{
    YAML::Node node = getNode(key);
    return node.as<string>();
}





/**
 * \brief      Return the absolute filename for the given parameter.
 * 
 * \details    The filename is first checked for a ENV['var'] pattern, which is then replaced by the 
 *             value of the environment value var.
 *             
 *             When the parameter contains an absolute filename, that value is returned.
 *             
 *             If the parameter contains a relative path, the filename is preceeded by the value of the 
 *             General/ProjectLocation parameter. When the ProjectLocation contains an environment 
 *             variable pattern ENV['var'], this is replaced with the value of the environment variable.
 *             
 *             An absolute path starts with a '/' character, otherwise the path is considered relative.
 *
 *             The key is the name of the input parameter. If the input parameter is part
 *             of a section or group, then the key is a combination of the group name and 
 *             the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 *             is part of the group ObservingParameters, then the key to get the value for this
 *             parameter would be "ObservingParameters/ExposureTime".
 * 
 * \param[in]  key The name of a parameter used in the PLATO Simulator
 *
 * \returns    An absolute filename
 */
string ConfigurationParameters::getAbsoluteFilename(const string &key) 
{
    using StringUtilities::replaceEnvironmentVariable;

    YAML::Node node = getNode(key);

	string filename = node.as<string>();

    filename = replaceEnvironmentVariable(filename);

	if (FileUtilities::isRelative(filename))
	{
		string projectLocation = this->getString("General/ProjectLocation");

        projectLocation = replaceEnvironmentVariable(projectLocation);

		return projectLocation + "/" + filename;
	}

	return filename;
}








/**
 * \brief      Set a (new) value to the given parameter.
 *
 * The parameter must exist, in which case it will be 
 * created automatically and assigned the given value.
 * 
 * The parameter might contain the name of a group/section in which it is defined. 
 * If this group does not exist, it will be created.
 * 
 * A Log message will be issued as a warning if the action overwrites an existing field.
 * 
 * \param[in]  key    The name of the field to assign the new value
 * \param[in]  value  The value to assign to the given field
 */
void ConfigurationParameters::setParameter(const string &key, const string &value)
{
    vector<string> fields = StringUtilities::split(key, '/');
    if (fields.size() > 1)
    {
        YAML::Node node = config[fields[0]];
        if (node)
        {
            YAML::Node subnode = node[fields[1]];
            if (subnode)
                Log.warning("ConfigurationParameters: setParameter() overwrites subnode \"" + fields[1] + "\" of node \"" + fields[0]);

            node[fields[1]] = value;
        }
        else
        {
            config[fields[0]][fields[1]] = value;
        }
    }
    else
    {
        if (config[key])
            Log.warning("ConfigurationParameters: setParameter() overwrites node \"" + key);
            
        config[key] = value;
    }
}






/**
 * \brief      Get the YAML node for the given key
 *
 * \details    
 * 
 * The key is the name of the input parameter. If the input parameter is part
 * of a section or group, then the key is a combination of the group name and 
 * the parameter name, separated by a '/' delimiter. E.g. if the parameter ExposureTime 
 * is part of the group ObservingParameters, then the key to get the value for this
 * parameter would be "ObservingParameters/ExposureTime". 
 * 
 * There can be multiple levels, e.g. "Camera/Distortion/Polynomial/Coefficients"
 * 
 * \param[in]  key   The name of a parameter used in the PLATO Simulator
 *
 * \return     The YAML node for the given key
 */
YAML::Node ConfigurationParameters::getNode(const string & key)
{
    vector<string> fields = StringUtilities::split(key, '/');

    // Why do we use a stack here?
    // Since we need to travers through the nodes and nodes are returned as references, so we can not
    // assign the returned node to the variable containing the previous node.
    // FIXME: explain this better, discuss with Joris...

    stack<YAML::Node> nodes;


    string parentNodeName = fields[0];

    nodes.push(config[parentNodeName]);

    if ( ! nodes.top() ) {
        string msg = "ConfigurationParameters: The field \"" + parentNodeName + "\" is not available in the configuration file (" + filename + ").";
        throw IllegalArgumentException(msg);
    }

    for (unsigned int idx = 1; idx < fields.size(); idx++)
    {
        nodes.push(nodes.top()[fields[idx]]);
        if ( ! nodes.top() ) {
            string msg = "ConfigurationParameters: The sub-field \"" + fields[idx] + "\" of field \"" + parentNodeName + "\" is not available in the configuration file (" + filename + ").";
            throw IllegalArgumentException(msg);
        }
        parentNodeName += "/" + fields[idx];
    }

    return nodes.top();
}




