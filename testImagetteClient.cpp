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
    zmq::socket_t socket (context, ZMQ_ROUTER);
    socket.bind ("tcp://*:5050");

    // initialize variables
    bool endOfSimulation = false;

    int imagetteCounter = 0;

    // this vector contains the identification of the client ,the imagette counter and whether the simulation is done
    std::vector<std::tuple<std::string, int, bool> >identityVec;

    bool lastImagetteFromLastClient = false;

    while(!lastImagetteFromLastClient)
    {
        // receive message from platosim
        zmq::message_t imagette;
        socket.recv (&imagette);

        std::string identity = std::string(static_cast<char*>(imagette.data()), imagette.size());

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
        socket.recv(&imagette);

        // get response from simulation
        socket.recv(&imagette);

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
            std::get<2>(*it) = true;
        }

        // get the rows and cols from the imagette message 
        int rows = receivedImagetteMessage.at(1);
        int cols = receivedImagetteMessage.at(2);

        // copy the vector without the first three values in an imagette vector
        std::vector<float> newImagette (receivedImagetteMessage.begin() + 3, receivedImagetteMessage.end());

        // print the imagette number
        std::cout << std::get<1>(*it) << std::endl;

/*        // print out the imagette
        for (int i = 0; i < cols; i++)
        {
            for (int j = 0; j < rows; j++)
            {
                std::cout << newImagette.at(i * (cols) + j) << "  ";
            }

            std::cout << std::endl;  
        }
*/

        // send identity

        zmq::message_t identityMessage(identity.size());
        memcpy (identityMessage.data(), identity.data(), identity.size());

        socket.send(identityMessage, ZMQ_SNDMORE);

        // send delimiter
        std::string delimiter = "";
        zmq::message_t delimiterMessage(delimiter.size());
        memcpy (delimiterMessage.data(), delimiter.data(), delimiter.size());

        socket.send(delimiterMessage, ZMQ_SNDMORE);



        // send answer back to platosim (this is not needed for platosim, but it is how zeroMQ works)
        std::string message = "imagette received";

        

        int strLength = message.length();

        //  Send reply back to client
        zmq::message_t reply (strLength);

        const char *cMessage = message.c_str(); 
        memcpy (reply.data (), cMessage, strLength);
        socket.send (reply);

        // end the sending of data when all clients have gotten the maximum number of steps
        auto boolIt = std::find_if(identityVec.begin(), identityVec.end(), [](const std::tuple<std::string, int, bool>& e) {return std::get<2>(e) == false;});

        if (boolIt == identityVec.end())
        {
            lastImagetteFromLastClient = true;
        }

    }



    return 0;
}