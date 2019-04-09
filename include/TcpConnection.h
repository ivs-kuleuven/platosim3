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
		TcpConnection(ConfigurationParameters &configParam, JitterGenerator* jitterFromNetwork, std::condition_variable* cond_var, std::mutex* m, bool* notified, bool* newStep);
		~TcpConnection();

		void connectToServer();
		void connectToClient();

	protected:

	private:

		void configure(ConfigurationParameters &configParams);
		std::vector <double> processServerReply(string replyString);

		bool endOfSimulation;
		string tcpAddress;

		JitterGenerator* jitterInstance;

		double internalTime;

		std::condition_variable* condVarPointer;
		std::mutex* mutexPointer;
		bool* notifiedPointer;
		bool* newStepPointer;
};