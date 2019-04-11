//
//  Binds REP socket to tcp://*:5555
//
#include <sstream>
#include <zmq.hpp>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#ifndef _WIN32
#include <unistd.h>
#else
#include <windows.h>


#define sleep(n)    Sleep(n)
#endif

int main () 
{
	//prepare random seed
	srand(time(NULL));

    //  Prepare our context and socket
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REP);
    socket.bind ("tcp://*:5555");

    int endOfSimulation = 0;

    double jitterStep = 0.0;
    double yaw;
    double pitch;
    double roll;

    int counter = 0;

    while (!endOfSimulation) 
    {
        if (counter == 100)
        {
            endOfSimulation = 1;
        }

        zmq::message_t request;

        //  Wait for next request from client
        socket.recv (&request);
        std::cout << "Received Hello" << std::endl;

        //  Do some 'work'
        sleep(1);

	// create some random jitter values

	yaw = (rand() % 100) / 10000.0;
	pitch = (rand() % 100) / 10000.0;
	roll = (rand() % 100) / 10000.0;


        std::ostringstream strEos; //, strJs, strY, strP, strR;
        std::ostringstream strJs;
        std::ostringstream strY;
        std::ostringstream strP;
        std::ostringstream strR;

        jitterStep += 0.1;

        strEos << endOfSimulation;
        strJs << jitterStep;
        strY << yaw;
        strP << pitch;
        strR << roll;

        std::string message = strEos.str() + "," + strJs.str() + "," + strY.str() + "," + strP.str() + "," + strR.str();

        std::cout << message << std::endl;

        int strLength = message.length();

        //  Send reply back to client
        zmq::message_t reply (strLength);

        const char *cMessage = message.c_str(); 

        memcpy (reply.data (), cMessage, strLength);
        socket.send (reply);

        counter++;
    }

    return 0;
}