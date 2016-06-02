#ifndef CONFIGURATION_PARAMETERS_H_
#define CONFIGURATION_PARAMETERS_H_

#include <string>
#include <regex>
#include <list>
#include <stack>
#include <fstream>

#include "yaml-cpp/yaml.h"
#include "yaml-cpp/exceptions.h"

#include "FileUtilities.h"
#include "StringUtilities.h"
#include "Exceptions.h"
#include "Logger.h"


using namespace std;


class ConfigurationParameters
{
    public:
        ConfigurationParameters();
        ConfigurationParameters(const string &);
        ~ConfigurationParameters();
    
        bool getBoolean(const string &);
        int getInteger(const string &);
        long getLong(const string &);
        double getDouble(const string &);
        string getString(const string &);
        string getAbsoluteFilename(const string &);
    
        vector <double> getDoubleVector(const string &key);
    
        void setParameter(const string &, const string &);
    
    private:
        string filename;
        YAML::Node config;
    
        YAML::Node getNode(const string & key);

};


#endif /* CONFIGURATION_PARAMETERS_H_ */
