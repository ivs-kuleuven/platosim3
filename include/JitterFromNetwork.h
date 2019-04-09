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

class JitterFromNetwork:public JitterGenerator
{
	public:

		JitterFromNetwork(ConfigurationParameters &configurationParameters, double readoutTimeBeforeNextExposure, std::condition_variable* cond_var, std::mutex* m, bool* notified, bool* newStep);

        	~JitterFromNetwork();

		virtual void configure(ConfigurationParameters &configParams);
		virtual tuple<double, double, double> getNextYawPitchRoll(double time) override;
        	virtual double getHeartbeatInterval() override;

		virtual void setCurrentJitterStep(double endSimulation, double timeStep, double yaw, double pitch, double roll) override;


	protected:


	private:

		std::condition_variable* condVarPointer;
		std::mutex* mutexPointer;
		bool* notifiedPointer;
		bool* newStepPointer;

		double oldTimeStep;
		double oldYaw;
		double oldPitch;
		double oldRoll;

		double currentTimeStep;
		double currentYaw;
		double currentPitch;
		double currentRoll;

};

#endif

