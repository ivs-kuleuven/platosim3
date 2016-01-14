#ifndef CONFIGURATION_PARAMETERS_H_
#define CONFIGURATION_PARAMETERS_H_

#include <string>
#include <list>

#include "ConfigurationFormat.h"


class ConfigurationParameters
{
public:
    ConfigurationParameters();
    ConfigurationParameters(const char*);
    ConfigurationParameters(const std::string &);
    ~ConfigurationParameters();

    // ConfigurationParameters getGroup(const std::string &);

    int getInteger(const std::string &);
    double getDouble(const std::string &);
    std::string getString(const std::string &);

private:
    ConfigurationFormat *format;
    std::string filename;
};


#endif /* CONFIGURATION_PARAMETERS_H_ */
