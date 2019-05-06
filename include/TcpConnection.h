#include "zmq.hpp"
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string>

#include "ConfigurationParameters.h"
#include "JitterGenerator.h"
#include "Detector.h"

class TcpConnection
{
	public:
		TcpConnection(ConfigurationParameters &configParam, std::condition_variable* cond_var, std::mutex* m, bool* notified, bool* newStep, bool* endSimulation);
		~TcpConnection();

		void connectToServer();
		void connectToClient();

		void setDetectorInstance(Detector* detector){detectorInstance = detector;};
		void setJitterInstance(JitterGenerator* jitter){jitterInstance = jitter;};

	protected:

	private:

		void configure(ConfigurationParameters &configParams);
		std::vector <double> processServerReply(string replyString);
		const char* convertMatrixToChar(arma::Mat<float>* pixelMapPointer);

		bool* endOfSimulation;
		string tcpAddressServer;
		string tcpAddressClient;

		JitterGenerator* jitterInstance;
		Detector* detectorInstance;

		double internalTime;

		std::condition_variable* condVarPointer;
		std::mutex* mutexPointer;
		bool* notifiedPointer;
		bool* newStepPointer;

		int numExposures;
};