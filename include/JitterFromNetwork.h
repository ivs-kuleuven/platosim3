#ifndef JITTERFROMNETWORK_H
#define JITTERFROMNETWORK_H

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



class JitterFromNetwork : public JitterGenerator
{
    public:

        static JitterGenerator* Instance(ConfigurationParameters &configParams);
        
        ~JitterFromNetwork();

        virtual void configure(ConfigurationParameters &configParams);
        virtual tuple<double, double, double> getNextYawPitchRoll(double time) override;
        virtual double getHeartbeatInterval() override;

    protected:

        JitterFromNetwork(ConfigurationParameters &configurationParameters);

        double yawRMS;              // [rad] 
        double pitchRMS;            // [rad]
        double rollRMS;             // [rad]
        double jitterTimeScale;     // [s]

        double lastYaw;
        double lastPitch;
        double lastRoll;

        double jitterTimeInterval;  // [s]

        long jitterNoiseSeed;

        double internalTime;        // [s]

    private:

        mt19937 jitterNoiseGenerator;
        normal_distribution<double> normalDistribution;

 };



#endif
 