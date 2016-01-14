#ifndef YAML_FORMAT_H_
#define YAML_FORMAT_H_

#include <string>

#include "ConfigurationFormat.h"
#include "yaml-cpp/yaml.h"

class YamlFormat : public ConfigurationFormat
{
public:
    YamlFormat(const std::string &);
    ~YamlFormat();

    int getInteger(const std::string &);
    double getDouble(const std::string &);
    std::string getString(const std::string &);

private:
    YAML::Node config;
};


#endif /* YAML_FORMAT_H_ */
