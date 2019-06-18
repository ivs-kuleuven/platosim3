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


std::string jitterFileName = std::getenv("PLATO_PROJECT_HOME") + std::string("/inputfiles/ohb18jun.txt");
int numberOfSteps = 80;
int numberOfSimulations = 2;


int main () 
{
    //  Prepare our context and socket
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_ROUTER);
    socket.bind ("tcp://*:5555");

    //create an ifstream
    std::ifstream infile(jitterFileName);

    double step, yaw, pitch, roll;
    int stepCounter = 0;

    std::vector<double> stepVec;
    std::vector<double> yawVec;
    std::vector<double> pitchVec;
    std::vector<double> rollVec;

    int endOfSimulation = 0;

    // read the file line by line and save the values in the respective vector
    while (infile >> step >> yaw >> pitch >> roll)
    {
        // make sure not to read more lines than used in the 
        if (stepCounter >= numberOfSteps)
        {
            break;
        }
        else
        {
            stepVec.push_back(step);
            yawVec.push_back(yaw);
            pitchVec.push_back(pitch);
            rollVec.push_back(roll);
            stepCounter++;
        }
    }

    // get the time between steps (assuming it won't change) in milliseconds
    
    // default, if there is only one entry in the file
    int stepLength = 1;

    if (stepVec.size() > 1)
    {
        stepLength = (stepVec.at(1) - stepVec.at(0)) * 1000;
    }


    bool lastStepToLastClient = false;

    // this vector contains the identification of the client ,the step counter and whether the simulation is done
    std::vector<std::tuple<std::string, int, bool> >identityVec;


    while(!lastStepToLastClient)
    {
        zmq::message_t request;


        // get identity
        socket.recv(&request);

        std::string identity = std::string(static_cast<char*>(request.data()), request.size());

        // search for identity in identity vector
        auto it = std::find_if(identityVec.begin(), identityVec.end(), [identity](const std::tuple<std::string, int, bool>& e) {return std::get<0>(e) == identity;});

        // if the identy is unknown create a new entry, if it is known, increment the counter
        if (it == identityVec.end())
        {
            identityVec.push_back(make_tuple(identity, 0, false));
            it = (identityVec.end() - 1);
        }
        else
        {
            std::get<1>(*it)++; 
        }
        
        // get delimiter
        socket.recv(&request);

        // get response from simulation
        socket.recv(&request);

        if (std::get<1>(*it) == stepVec.size()-1)
        {
            std::cout << "end of simulation" << std::endl; 
            std::get<2>(*it) = true;
        }

        // create the next message
        std::ostringstream strEos;
        std::ostringstream strJs;
        std::ostringstream strY;
        std::ostringstream strP;
        std::ostringstream strR;

        strEos << std::get<2>(*it);
        strJs << stepVec.at(std::get<1>(*it));
        strY << yawVec.at(std::get<1>(*it));
        strP << pitchVec.at(std::get<1>(*it));
        strR << rollVec.at(std::get<1>(*it));

        std::string message = strEos.str() + "," + strJs.str() + "," + strY.str() + "," + strP.str() + "," + strR.str();

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

        // send jitterStep
        zmq::message_t reply (strLength);

        const char *cMessage = message.c_str(); 

        memcpy (reply.data (), cMessage, strLength);
        socket.send (reply);
        std::cout << "Sent next jitter step to: " << identity << std::endl;

        // end the sending of data when all clients have gotten the maximum number of steps
        auto boolIt = std::find_if(identityVec.begin(), identityVec.end(), [](const std::tuple<std::string, int, bool>& e) {return std::get<2>(e) == false;});

        if (boolIt == identityVec.end())
        {
            lastStepToLastClient = true;
        }

    }

    return 0;
}