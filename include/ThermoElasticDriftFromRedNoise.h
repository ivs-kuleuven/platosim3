#ifndef THERMOELASTICDRIFTFROMREDNOISE_H
#define THERMOELASTICDRIFTFROMREDNOISE_H

#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>

#include "Logger.h"
#include "Units.h"
#include "Exceptions.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "DriftGenerator.h"


using namespace std;



class ThermoElasticDriftFromRedNoise : public DriftGenerator
{
    public:

        ThermoElasticDriftFromRedNoise(ConfigurationParameters &configurationParameters, double readoutTimeBeforeNextExposure);
        ~ThermoElasticDriftFromRedNoise();

        virtual void configure(ConfigurationParameters &configParams, double readoutTimeBeforeNextExposure);
        virtual tuple<double, double, double> getNextYawPitchRoll(double time) override;
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

        double internalTime;        // [s]

    private:

        mt19937 driftNoiseGenerator;
        normal_distribution<double> normalDistribution;

 };



#endif
 
