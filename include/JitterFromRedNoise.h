#ifndef JITTERFROMREDNOISE_H
#define JITTERFROMREDNOISE_H

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
#include "JitterGenerator.h"


using namespace std;



class JitterFromRedNoise : public JitterGenerator
{
    public:

        JitterFromRedNoise(ConfigurationParameters &configurationParameters);
        ~JitterFromRedNoise();

        virtual void configure(ConfigurationParameters &configParams);
        virtual tuple<double, double, double> getNextYawPitchRoll(double timeInterval) override;
        virtual double getHeartbeatInterval() override;

    protected:

        double yawRMS;              // [rad] 
        double pitchRMS;            // [rad]
        double rollRMS;             // [rad]
        double jitterTimeScale;     // [s]

        double lastYaw;
        double lastPitch;
        double lastRoll;

        double jitterTimeInterval;  // [s]

        long jitterNoiseSeed;

    private:

        mt19937 jitterNoiseGenerator;
        normal_distribution<double> normalDistribution;

 };



#endif
 