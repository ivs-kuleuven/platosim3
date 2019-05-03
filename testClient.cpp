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


#endif

int main () 
{
    //  Prepare our context and socket
    zmq::context_t context (1);
    zmq::socket_t socket (context, ZMQ_REP);
    socket.bind ("tcp://*:5050");

    // initialize variables
    bool endOfSimulation = false;

    int imagetteCounter = 0;

    // repeat until the simulation is at an end
    while (!endOfSimulation) 
    {
        // receive message from platosim
        zmq::message_t imagette;
        socket.recv (&imagette);
        std::cout << "Received Notification from PlatoSim" << std::endl;

        // convert the message from platosim to string
        std::string sImagette = std::string(static_cast<char*>(imagette.data()), imagette.size());
        
        // initiate a stringstream
        std::stringstream ss(sImagette);

        // initiate a vector of doubles to write the the imagette values into
        std::vector<double> receivedImagetteMessage;

        double i;

        // convert the string to double values

        while (ss >> i)
        {
            receivedImagetteMessage.push_back(i);

            if (ss.peek() == ',' || ss.peek() == ' ')
            {
                ss.ignore();
            }
        }

        // check whether the simulation is to end
        if (receivedImagetteMessage.at(0) == 1)
        {
            endOfSimulation = true;
        }

        // get the rows and cols from the imagette message 
        int rows = receivedImagetteMessage.at(1);
        int cols = receivedImagetteMessage.at(2);

        // copy the vector without the first three values in an imagette vector
        std::vector<float> newImagette (receivedImagetteMessage.begin() + 3, receivedImagetteMessage.end());

        // print the imagette number
        std::cout << imagetteCounter << std::endl;

        // print out the imagette
        for (int i = 0; i < cols; i++)
        {
            for (int j = 0; j < rows; j++)
            {
                std::cout << newImagette.at(i * (cols) + j) << "  ";
            }

            std::cout << std::endl;   
        }

        // increase the counter
        imagetteCounter++;

        // send answer back to platosim (this is not needed for platosim, but it is how zeroMQ works)
        std::string message = "imagette received";

        int strLength = message.length();

        //  Send reply back to client
        zmq::message_t reply (strLength);

        const char *cMessage = message.c_str(); 
        memcpy (reply.data (), cMessage, strLength);
        socket.send (reply);
    }

    return 0;
}