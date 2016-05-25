#ifndef THERMOELASTICDRIFTFROMREDNOISE_H
#define THERMOELASTICDRIFTFROMREDNOISE_H

#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>

#include "Logger.h"
#include "Units.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "DriftGenerator.h"


using namespace std;



class ThermoElasticDriftFromRedNoise : public DriftGenerator
{
    public:

        ThermoElasticDriftFromRedNoise(ConfigurationParameters &configurationParameters);
        ~ThermoElasticDriftFromRedNoise();

        virtual void configure(ConfigurationParameters &configParams);
        virtual tuple<double, double, double> getNextYawPitchRoll(double timeInterval) override;
        virtual double getHeartbeatInterval() override;

    protected:

        double yawRMS;              // [rad] 
        double pitchRMS;            // [rad]
        double rollRMS;             // [rad]
        double driftTimeScale;      // [s]

        double lastYaw;
        double lastPitch;
        double lastRoll;

        double driftTimeInterval;   // [s]

        long driftNoiseSeed;

    private:

        mt19937 driftNoiseGenerator;
        normal_distribution<double> normalDistribution;

 };



#endif
 