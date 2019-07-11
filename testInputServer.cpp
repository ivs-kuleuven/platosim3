//
//  Binds REP socket to tcp://*:5555
//

#include <chrono>
#include <thread>
#include <sstream>
#include <fstream>
#include <zmq.hpp>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <chrono>
#include <algorithm>

#ifndef _WIN32
#include <unistd.h>
#else
#include <windows.h>

#endif


// set the path to a file with the subfield positions and size for the uses stars
std::string serverInputFileName = std::getenv("PLATO_PROJECT_HOME") + std::string("/starPositions.txt");

// initialize after how many exposures the subfield size will change
int changeSizeAtExposure = 20000;

int main () 
{
    // create a data vector

    std::vector<std::tuple<int, int, int, int, int, int, int, int, int, double, double, int>> starPositionVec;

    int starId, subFieldRows1, subFieldCols1, subFieldRows2, subFieldCols2, rowPos1, rowPos2, colPos1, colPos2, orientation;

    double offsetX, offsetY;

    //create an ifstream
    std::ifstream infile(serverInputFileName); 

    // read all information from the file to a vector
    while (infile >> starId >> subFieldCols1 >> subFieldRows1 >> colPos1 >> rowPos1 >> subFieldCols2 >> subFieldRows2 >> colPos2 >> rowPos2 >> offsetX >> offsetY >> orientation)
    {
        starPositionVec.push_back(std::make_tuple(starId, subFieldRows1, subFieldCols1, rowPos1, colPos1, subFieldRows2, subFieldCols2, rowPos2, colPos2, offsetX, offsetY, orientation));
    }


    //  Prepare our context and socket
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_ROUTER);
    socket.bind ("tcp://*:5000");

    //int numSimulation = 0;

    bool lastStepToLastClient = false;

    // this vector contains the identification of the client, the exposure counter and whether the simulation is done
    std::vector<std::tuple<std::string, int, bool> >identityVec;

    while(!lastStepToLastClient)
    {
        zmq::message_t platoMessage;

        // get identity from the message
        socket.recv(&platoMessage);

        std::string identity = std::string(static_cast<char*>(platoMessage.data()), platoMessage.size());
        
        // get delimiter
        socket.recv(&platoMessage);

        // get response from simulation
        socket.recv(&platoMessage);

        // search for identity in identity vector
        auto it = std::find_if(identityVec.begin(), identityVec.end(), [ = ](auto item) {return std::get<0>( item ) == identity;});

        // if the identy is unknown create a new entry and assign it a vector iterator
        if (it == identityVec.end())
        {
            identityVec.push_back(std::make_tuple(identity, 0, false));

            std::cout << "New simulation message detetcted. Simulation Nr: " << identityVec.size() << std::endl; 

            it = (identityVec.end() - 1);
        }

        // convert the message from platosim to string
        std::string sPlatoMessage = std::string(static_cast<char*>(platoMessage.data()), platoMessage.size());
        
        // initiate a stringstream
        std::stringstream ss(sPlatoMessage);

        // initiate a vector of doubles to write the the imagette values into
        std::vector<double> receivedPlatoMessage;

        double i;

        // convert the string to double values

        while (ss >> i)
        {
            receivedPlatoMessage.push_back(i);

            if (ss.peek() == ',' || ss.peek() == ' ')
            {
                ss.ignore();
            }
        }

        // set the identity vector with exposureNumber and endOfSimulation
        std::get<1>(*it) = receivedPlatoMessage.at(1);
        std::get<2>(*it) = receivedPlatoMessage.at(0);

        int strSubFieldRows;
        int strSubFieldCols;
        int strRowPos;
        int strColPos;
        double strOffsetX;
        double strOffsetY;
        int strOrientation;

        // print the number of exposures and whether the simulation is at an end
        std::cout << std::to_string(std::get<1>(*it)) << std::endl;
        std::cout << std::to_string(std::get<2>(*it)) << std::endl;

        // dependend on how many exposures are already done, take respective size and position of the vector entry
        if (std::get<1>(*it) < changeSizeAtExposure)
        {
            strSubFieldCols = std::get<1>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strSubFieldRows = std::get<2>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strColPos = std::get<3>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strRowPos = std::get<4>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strOffsetX = std::get<9>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strOffsetY = std::get<10>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strOrientation = std::get<11>(starPositionVec.at(std::distance(identityVec.begin(), it)));
        }
        else
        {
            strSubFieldCols = std::get<5>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strSubFieldRows = std::get<6>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strColPos = std::get<7>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strRowPos = std::get<8>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strOffsetX = std::get<9>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strOffsetY = std::get<10>(starPositionVec.at(std::distance(identityVec.begin(), it)));
            strOrientation = std::get<11>(starPositionVec.at(std::distance(identityVec.begin(), it)));
        }

        // compose a message to the simulation client
        std::string message = std::to_string(strSubFieldRows) + "," + std::to_string(strSubFieldCols) + "," + std::to_string(strColPos) + "," + std::to_string(strRowPos) + "," + std::to_string(strOffsetX) + "," + std::to_string(strOffsetY) + "," + std::to_string(strOrientation);

        std::cout << message << std::endl;

        int strLength = message.length();

        // send identity

        zmq::message_t identityMessage(identity.size());
        memcpy (identityMessage.data(), identity.data(), identity.size());

        socket.send(identityMessage, ZMQ_SNDMORE);

        // send delimiter
        std::string delimiter = "";
        zmq::message_t delimiterMessage(delimiter.size());
        memcpy (delimiterMessage.data(), delimiter.data(), delimiter.size());

        socket.send(delimiterMessage, ZMQ_SNDMORE);

        // send input variables
        zmq::message_t reply (strLength);

        const char *cMessage = message.c_str(); 

        memcpy (reply.data (), cMessage, strLength);
        socket.send (reply);
        std::cout << "Sent input parameters to: " << identity << std::endl;

        // end the sending of data when all clients have gotten the maximum number of steps
        auto boolIt = std::find_if(identityVec.begin(), identityVec.end(), [ = ](auto item) {return std::get<2>(item) == false;});

        if (boolIt == identityVec.end())
        {
            lastStepToLastClient = true;
        }

    }

    return 0;
}