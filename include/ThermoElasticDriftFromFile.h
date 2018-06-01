#ifndef THERMOELASTICDRIFTFROMFILE_H
#define THERMOELASTICDRIFTFROMFILE_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

#include "Exceptions.h"
#include "Logger.h"
#include "Units.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "DriftGenerator.h"


using namespace std;



class ThermoElasticDriftFromFile : public DriftGenerator
{
    public:

        static DriftGenerator* Instance(ConfigurationParameters &configParams);

        ~ThermoElasticDriftFromFile();

        virtual void configure(ConfigurationParameters &configParams);
        virtual tuple<double, double, double> getNextYawPitchRoll(double time) override;
        virtual double getHeartbeatInterval() override;

    protected:

        ThermoElasticDriftFromFile(ConfigurationParameters &configurationParameters);

        string pathToDriftFile;
        double beginTime;                 // Only read the drift file from beginTime to endTime
        double endTime;

        int timeIndex;
        double internalTime;

        vector<double> timeFromFile;    // [s]
        vector<double> yaw;             // [rad]
        vector<double> pitch;           // [rad]
        vector<double> roll;            // [rad]


    private:

 };



#endif
 
