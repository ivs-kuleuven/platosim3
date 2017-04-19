#ifndef INCLUDE_TEMPERATUREFROMFILE_H_
#define INCLUDE_TEMPERATUREFROMFILE_H_

#include <string>
#include <vector>

#include "Logger.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "TemperatureGenerator.h"

class TemperatureFromFile : public TemperatureGenerator
{
    public:

		TemperatureFromFile(ConfigurationParameters &configurationParameters, string component);
        ~TemperatureFromFile();

        virtual void configure(ConfigurationParameters &configParams, string component);
        virtual double getNextTemperature(double time);

    protected:

        string pathToTemperatureFile;

        int timeIndex;

        vector<double> timeFromFile;    // [s]
        vector<double> temperature;     // [K]


    private:

 };




#endif /* INCLUDE_TEMPERATUREFROMFILE_H_ */
