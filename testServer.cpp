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

#ifndef _WIN32
#include <unistd.h>
#else
#include <windows.h>

#endif


std::string jitterFileName = std::getenv("PLATO_PROJECT_HOME") + std::string("/inputfiles/ohb18jun.txt");
int numberOfSteps = 80;

int main () 
{
    //  Prepare our context and socket
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REP);
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

    std::vector<double> timeVec;
    double completeTime = 0;

    for (int i = 0; i < numberOfSteps; i++)
    {
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        zmq::message_t request;

        //  Wait for next request from client
        socket.recv (&request);
        std::cout << "Received request for next jitter step" << std::endl;

        // Do some 'work'
        std::this_thread::sleep_for(std::chrono::milliseconds(stepLength));

	if (i == stepVec.size()-1)
        {
            std::cout << "end of simulation" << std::endl; 
            endOfSimulation = 1;
        }

        std::ostringstream strEos;
        std::ostringstream strJs;
        std::ostringstream strY;
        std::ostringstream strP;
        std::ostringstream strR;

        strEos << endOfSimulation;
        strJs << stepVec.at(i);
        strY << yawVec.at(i);
        strP << pitchVec.at(i);
        strR << rollVec.at(i);

        std::string message = strEos.str() + "," + strJs.str() + "," + strY.str() + "," + strP.str() + "," + strR.str();

        std::cout << message << std::endl;

        int strLength = message.length();

        //  Send reply back to client
        zmq::message_t reply (strLength);

        const char *cMessage = message.c_str(); 

        memcpy (reply.data (), cMessage, strLength);
        socket.send (reply);

	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
        
        double timeDiff = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() /1000000.0;

	timeVec.push_back(timeDiff);

        completeTime += timeDiff;

    }

    

    std::cout << *max_element(std::begin(timeVec), std::end(timeVec)) << std::endl;

    std::cout << completeTime/timeVec.size() << std::endl;

    return 0;
}