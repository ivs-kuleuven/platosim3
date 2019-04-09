#include "zmq.hpp"
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string>

#include "ConfigurationParameters.h"
#include "JitterGenerator.h"

class TcpConnection
{
	public:
		TcpConnection(ConfigurationParameters &configParam, JitterGenerator* jitterFromNetwork);
		~TcpConnection();

		void connectToServer();
		void connectToClient();

	protected:

	private:

		void configure(ConfigurationParameters &configParams);

		bool endOfSimulation;

		JitterGenerator* jitterInstance;
};