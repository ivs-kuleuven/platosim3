#ifndef INCLUDE_TEMPERATUREGENERATOR_H_
#define INCLUDE_TEMPERATUREGENERATOR_H_


#include <string>
#include <vector>

#include "Logger.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"

using namespace std;



class TemperatureGenerator
{
    public:

		TemperatureGenerator(){};
        virtual ~TemperatureGenerator(){};

        virtual double getNextTemperature(double time) = 0;

    protected:

        double internalTime = 0.0;

    private:

 };




#endif /* INCLUDE_TEMPERATUREGENERATOR_H_ */
