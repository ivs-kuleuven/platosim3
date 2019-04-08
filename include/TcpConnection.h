#include "zmq.hpp"
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string>

#include "ConfigurationParameters.h"


class TcpConnection
{
	public:
		TcpConnection(string inputFileName);
		~TcpConnection();

	protected:

	private:

		void configure(ConfigurationParameters &configParams);
}