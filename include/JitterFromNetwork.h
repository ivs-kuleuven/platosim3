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

#include "zmq.hpp"


using namespace std;



class JitterFromNetwork : public JitterGenerator
{
    public:

        JitterFromNetwork(ConfigurationParameters &configurationParameters, zmq::socket_t* socketPtr);
        ~JitterFromNetwork();

        virtual void configure(ConfigurationParameters &configParams);
        virtual tuple<double, double, double> getNextYawPitchRoll(double time) override;
        virtual double getHeartbeatInterval() override;

        virtual std::tuple<bool, double, double, double, double> processJitterMessage(std::string message);

    protected:

        double currentTimeStep;     // [s]

        double currentYaw;             // [rad]
        double currentPitch;           // [rad]
        double currentRoll;            // [rad]

        double oldTimeStep;        // [s]

        double oldYaw;             // [rad]
        double oldPitch;           // [rad]
        double oldRoll;            // [rad]

        double jitterTimeInterval;  // [s]

        double internalTime;        // [s]

        zmq::socket_t* jitterSocketPtr;

    private:

 };



#endif
 
