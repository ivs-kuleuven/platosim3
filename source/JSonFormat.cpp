#include <string>
#include <fstream>

#include "StringUtilities.h"
#include "JSonFormat.h"
#include "json.hpp"

using json = nlohmann::json;

JSonFormat::JSonFormat(const std::string &name)
{
        std::ifstream ifs (name);
        ifs >> j_input;
}




int JSonFormat::getInteger(const std::string &key) 
{
    std::vector<std::string> fields;

    split(fields, key, "/", split::no_empties);

    if (fields.size() > 1)
    {
        json node = j_input[fields[0]];
        return node[fields[1]].get<int>();
    }
    else 
    {
        json j_value = j_input[key];
        return j_value.get<int>();
    }
}

double JSonFormat::getDouble(const std::string &key) 
{
    json j_value = j_input[key];
    return j_value.get<double>();
}

std::string JSonFormat::getString(const std::string &key) 
{
    json j_value = j_input[key];
    return j_value.get<std::string>();
}




JSonFormat::~JSonFormat() {}
