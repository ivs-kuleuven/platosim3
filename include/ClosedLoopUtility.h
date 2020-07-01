#ifndef CLOSEDLOOPUTILITY_H
#define CLOSEDLOOPUTILITY_H


#include <string>

#include "zmq.hpp"
#include "ConfigurationParameters.h"
#include "armadillo"
#include "Logger.h"





// this class governs the sending and receiving of messages between the used 
// detector class and the imagette client and window position server

class ClosedLoopUtility
{
    
    public:

        ClosedLoopUtility(ConfigurationParameters &configParams);

        ~ClosedLoopUtility() {};

    protected:

        std::tuple<bool, uint, uint, uint, uint, double> getNewWindowPosition(double exposureTime);

        void sendImagetteToClient(arma::Mat<float>* imagette, uint exposureNumber);

        bool firstExposure;
        bool sendImagettesToClient;
        bool getWindowPositionFromServer;


    private:

        virtual void configure(ConfigurationParameters &configParams);

        void sendMessageToSocket(const std::string messageString, zmq::socket_t* socketPointer);

        std::string receiveMessageFromSocket(zmq::socket_t* socketPointer);

        std::vector<double> convertStringToDoubleVec(std::string message);


        std::string convertMatrixToString(arma::Mat<float>* pixelMapPointer, const uint exposureNumber);

        int windowPositionSocketTimeout;

        string windowPositionAddress;
        string imagetteAddress;

        zmq::context_t context;

        zmq::socket_t imagetteSocket;
        zmq::socket_t windowPositionSocket;
};



#endif
