#ifndef INCLUDE_NOMINALTEMPERATURE_H_
#define INCLUDE_NOMINALTEMPERATURE_H_

#include <string>
#include <vector>

#include "Logger.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "TemperatureGenerator.h"

class NominalTemperature : public TemperatureGenerator
{
    public:

		NominalTemperature(ConfigurationParameters &configurationParameters, string component);
        ~NominalTemperature();

        virtual void configure(ConfigurationParameters &configParams, string component);
        virtual double getNextTemperature(double time);

    protected:

        double nominalTemperature;


    private:

 };


#endif /* INCLUDE_NOMINALTEMPERATURE_H_ */
