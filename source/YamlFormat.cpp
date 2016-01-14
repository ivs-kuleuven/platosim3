#include <string>

#include "YamlFormat.h"
#include "StringUtilities.h"
#include "Exceptions.h"

YamlFormat::YamlFormat(const std::string &name)
{
    config = YAML::LoadFile(name);
}



int YamlFormat::getInteger(const std::string &key)
{
    std::vector<std::string> fields;

    split(fields, key, "/", split::no_empties);

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

double YamlFormat::getDouble(const std::string &key) 
{
    return config[key].as<double>();
}

std::string YamlFormat::getString(const std::string &key) 
{
    return config[key].as<std::string>();
}


YamlFormat::~YamlFormat() {}
