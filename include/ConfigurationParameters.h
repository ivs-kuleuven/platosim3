#ifndef CONFIGURATION_PARAMETERS_H_
#define CONFIGURATION_PARAMETERS_H_

#include <string>
#include <list>
#include <fstream>

#include "yaml-cpp/yaml.h"

#include "FileUtilities.h"
#include "StringUtilities.h"
#include "Exceptions.h"
#include "Logger.h"


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
    long getLong(const string &);
    double getDouble(const string &);
    string getString(const string &);
    string getAbsoluteFilename(const string &);

    void setParameter(const string &, const string &);

private:
    string filename;
    YAML::Node config;
};


#endif /* CONFIGURATION_PARAMETERS_H_ */
