#ifndef CONFIGURATION_PARAMETERS_H_
#define CONFIGURATION_PARAMETERS_H_

#include <string>
#include <fstream>

#include "Logger.h"
#include "yaml-cpp/yaml.h"

using namespace std;


class ConfigurationParameters
{
public:
    ConfigurationParameters();
    ConfigurationParameters(const char*);
    ConfigurationParameters(const string &);
    ~ConfigurationParameters();

    bool getBoolean(const string &);
    int getInteger(const string &);
    double getDouble(const string &);
    string getString(const string &);
    string getAbsoluteFileName(const string &);

    void setNode(const string &, const string &);

private:
    string filename;
    YAML::Node config;
};


#endif /* CONFIGURATION_PARAMETERS_H_ */
