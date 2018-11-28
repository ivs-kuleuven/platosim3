#include "zmq.hpp"
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string>

#include "ConfigurationParameters.h"

#include "JitterGenerator.h"
#include "NoJitter.h"
#include "JitterFromFile.h"
#include "JitterFromRedNoise.h"
#include "JitterFromNetwork.h"

#include "Simulation.h"

class TcpConnection
{
	public:
		TcpConnection(string inputFilename);
		~TcpConnection();

		void connectToServer(Simulation* simulationInstance, std::condition_variable* cond_var, bool *notified, bool* newStep, std::mutex* m);
		void connectToClient();

	protected:

	private:

		void configure(ConfigurationParameters &configParams);

		bool newJitterStep;

		bool useJitter;
		string jitterSource;

		double internalTime;

		bool endOfSimulation;

		std::vector<double> currentJitterStepVec;

		bool client;
        bool server;

        std::string tcpAddress;
};