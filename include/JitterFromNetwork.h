#ifndef JITTERFROMNETWORK_H
#define JITTERFROMNETWORK_H

#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>
#include <condition_variable>

#include "zmq.hpp"
#include <chrono>
#include <thread>

#include "Logger.h"
#include "Units.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "JitterGenerator.h"



using namespace std;

class TcpConnection;

class JitterFromNetwork : public JitterGenerator
{
    public:

        static JitterGenerator* Instance(ConfigurationParameters &configurationParameters, double readoutTimeBeforeNextExposure);
        static JitterGenerator* Instance(ConfigurationParameters &configurationParameters, double readoutTimeBeforeNextExposure, TcpConnection* tcpConnection);
        
        ~JitterFromNetwork();

        virtual void configure(ConfigurationParameters &configParams);
        virtual tuple<double, double, double> getNextYawPitchRoll(double time) override;
        virtual double getHeartbeatInterval() override;

        virtual void setCurrentJitterStep(double endSimulation, double timeStep, double yaw, double pitch, double roll) override;

        virtual bool isClient() override;

        virtual bool simulationEnd() override;

    protected:

        JitterFromNetwork(ConfigurationParameters &configurationParameters, double readoutTimeBeforeNextExposure);
        JitterFromNetwork(ConfigurationParameters &configurationParameters, double readoutTimeBeforeNextExposure, TcpConnection* tcpConnection);

        double yawRMS;              // [rad] 
        double pitchRMS;            // [rad]
        double rollRMS;             // [rad]
        double jitterTimeScale;     // [s]

        double lastYaw;
        double lastPitch;
        double lastRoll;

        long jitterNoiseSeed;

        double internalTime;        // [s]

        string serverAddress;



    private:

        mt19937 jitterNoiseGenerator;
        normal_distribution<double> normalDistribution;

        TcpConnection* tcpConnectionInstance;

        double oldTimeStep;
        double oldYaw;
        double oldPitch;
        double oldRoll;

        double currentTimeStep;
        double currentYaw;
        double currentPitch;
        double currentRoll;

        double jitterTimeInterval;  // [s]

        bool endSimulation;
 };



#endif
 