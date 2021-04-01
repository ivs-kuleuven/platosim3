#include <zmq.hpp>
#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <thread>
#include <unistd.h>


typedef unsigned int uint;

bool endTest = false;
bool simsStarted = false;
bool simsRunning = false;
bool simsReady = false;
uint simsCreated = 0;
uint simsCheckedIn = 0;

std::string platoSimWinPositionAddress = "tcp://*:5558";
std::string platoSimJitterAddress = "tcp://*:5559";
std::string platoSimImagetteAddress = "tcp://*:5560";

std::string jitterFileName = std::string(std::getenv("PLATO_PROJECT_HOME")) + "/testData/shortJitter.txt";

std::string winPositionFileName = std::string(std::getenv("PLATO_PROJECT_HOME")) + "/testData/windowPositions.txt";


std::vector<std::tuple<int, int, int, int, int, int, int, int, int, int>> starPositionVec;


int stepCounter = 0;

uint simNum = 0;

std::vector<double> stepVec;
std::vector<double> yawVec;
std::vector<double> pitchVec;
std::vector<double> rollVec;

std::vector<std::string> identityVec;


void sendMessageToSocket(std::string messageString, std::string identityString, zmq::socket_t* socketPointer)
{
    if (identityString.length() != 0)
    {
        zmq::message_t identity(identityString.length());

        const char *cIdentity = identityString.c_str(); 

        memcpy (identity.data (), cIdentity, identityString.length());

        socketPointer->send(identity, ZMQ_SNDMORE);
    }

    zmq::message_t message(messageString.length());

    const char *cMessage = messageString.c_str(); 

    memcpy (message.data (), cMessage, messageString.length());

    socketPointer->send(message);
}



void receiveMessages(zmq::socket_t* imagetteSocket, zmq::socket_t* winPositionSocket)
{
    // initialize a poll set
    zmq::pollitem_t socketConnections[] =
    {
        {static_cast<void*> (*imagetteSocket), 0, ZMQ_POLLIN, 0},
        {static_cast<void*> (*winPositionSocket), 0, ZMQ_POLLIN, 0}
    };    

    while(!endTest)
    {
        std::cout << "Waiting for PlatoSim Messages" << std::endl;

        zmq::message_t message;

        zmq::poll (&socketConnections [0], 2, -1);

        // check wether there are new messages over the platosim imagette socket
        if(socketConnections[0].revents & ZMQ_POLLIN)
        {
            // get the identity of the message delivered
            
            imagetteSocket->recv(&message);

            std::string identityString = std::string(static_cast<char*>(message.data()), message.size());

            std::cout << identityString << std::endl;

            // get the imagette message

            imagetteSocket->recv(&message);

            std::string messageString = std::string(static_cast<char*>(message.data()), message.size());


            std::vector<std::string> receivedMessage;

            // convert the message string to single strings within a vector

            std::size_t prev = 0, pos;
            while ((pos = messageString.find_first_of(" ,;", prev)) != std::string::npos)
            {
                if (pos > prev)
                {
                    receivedMessage.push_back(messageString.substr(prev, pos - prev));
                }
                prev = pos+1;
            }
            if (prev < messageString.length())
            {
                receivedMessage.push_back(messageString.substr(prev, std::string::npos));
            }

            uint cols = std::stoi(receivedMessage.at(2));

            std::vector<std::string>::const_iterator first = receivedMessage.begin() + 3;
            std::vector<std::string>::const_iterator last = receivedMessage.end();
            std::vector<std::string> messageBody(first, last);

            std::string imagetteString;

            if (!endTest)
            {
                // convert and print the imagette
                for (uint i = 0; i < messageBody.size(); i += cols)
                {
                    for (uint j = 0; j < cols; j++)
                    {
                        std::cout << messageBody.at(i + j) << " ";

                        imagetteString += messageBody.at(i + j) + " ";
                    }

                    std::cout << std::endl;

                }
            }

        }


        // check wether there are new messages over the platosim imagette socket
        if(socketConnections[1].revents & ZMQ_POLLIN)
        {
            // get the identity of the message delivered
            
            winPositionSocket->recv(&message);

            std::string identityString = std::string(static_cast<char*>(message.data()), message.size());

            identityVec.push_back(identityString);

            // get the imagette message

            winPositionSocket->recv(&message);

            std::string messageString = std::string(static_cast<char*>(message.data()), message.size());

            std::cout << "Got message from Simulation Instance: " << identityString << std::endl;

            std::cout << "Message: " << messageString << std::endl;

            simsCheckedIn++;



            if (simsCheckedIn == simsCreated)
            {
                simsReady = true;

                std::cout << "All simualtions checked in." << std::endl;
 
            }
        }
    }
}


void sendMessages(zmq::socket_t* jitterSocket, zmq::socket_t* winPositionSocket)
{
    while(!endTest)
    {
        std::cout << std::endl;

        std::cout << "Which TC would you like to send?" << std::endl << std::endl;

        std::cout << "Start n Simulations:                  1" << std::endl;

        std::cout << "Send window positions to Simulations: 2" << std::endl;

        std::cout << "Send 20 jitter steps:                 3" << std::endl;

        std::cout << "Change window size:                   4" << std::endl;

        std::cout << "Stop Simulations:                     5" << std::endl;

        std::cout << "Remove output files:                  6" << std::endl;

        std::cout << std::endl;

        std::cout << "End Program:                          0" << std::endl << std::endl;


        char c;

        std::cin >> c;

        uint input;

        try
        {
            input = (int)c - 48;

        }
        catch (...)
        {
            std::cout << "Please enter a valid number between 0 and 5" << std::endl;

            continue;  
        }

        // end the programm
        if (input == 0)
        {
            // if no simulation is running create one (this is necessary to wake up the receive messages thread, so it can shut down)
            if (!simsRunning)
            {
                std::string systemCallString = "xterm -e ./platosim $PLATO_PROJECT_HOME/inputfiles/inputfgs.yaml out.hdf5 log.txt &";

                const char* systemCall = systemCallString.c_str(); 

                system(systemCall);

                usleep(1000000);

            }

            for (uint i = 0; i < identityVec.size(); i++)
            {
                sendMessageToSocket("6 6 0 0 0", identityVec.at(i), winPositionSocket);

                sendMessageToSocket("", identityVec.at(i), jitterSocket);
            }



            if (!simsRunning)
            {
                std::string removeString = "rm ./out.hdf5 ./log.txt";

                const char* removeCall = removeString.c_str(); 

                system(removeCall);
            }

            endTest = true;

        }
        // start n simulations
        else if (input == 1)
        {
            std::cout << "How many Simualtions du you want to start? (1-5)" << std::endl;

            char d;

            std::cin >> d;

            try
            {
                simNum = (int)d - 48;
            }
            catch (...)
            {
                std::cout << "Please enter a valid number between 1 and 5" << std::endl;

                continue;
            }

            if (simNum <= 0 || simNum > 6)
            {
                std::cout << "Please enter a valid number between 1 and 5" << std::endl;

                continue;
            }
            
            if (simsStarted)
            {
                std::cout << "Simulations are already running" << std::endl;

                continue;
            }
            else
            {

                for (uint i = 0; i < simNum; i++)
                {

                    std::string systemCallString = "xterm -e ./platosim $PLATO_PROJECT_HOME/inputfiles/inputfgs.yaml out" + std::to_string(i) + ".hdf5 log" + std::to_string(i) + ".txt &";

                    const char* systemCall = systemCallString.c_str(); 

                    system(systemCall);

                    simsCreated++;

                }

                simsStarted = true;
            }

        }
        // send win position message
        else if (input == 2)
        {
            if(simsReady)
            {
                std::cout << "Send the window positions to the simulations" << std::endl;

                for (uint i = 0; i < simNum; i++)
                {
                    std::cout << "Send the window position to identity: " << identityVec.at(i) << std::endl;
                
                    std::string identityString = identityVec.at(i);

                    std::tuple<int, int, int, int, int, int, int, int, int, int> winPos = starPositionVec.at(0);

                    std::string winPosString =  std::to_string(std::get<1>(winPos)) + " " +
                                                std::to_string(std::get<2>(winPos)) + " " +
                                                std::to_string(std::get<3>(winPos)) + " " +
                                                std::to_string(std::get<4>(winPos)) + " " +
                                                std::to_string(std::get<9>(winPos));

                    std::cout << "Identity: " << identityString << std::endl;

                    std::cout << "Message: " << winPosString << std::endl;

                    // send position string to identity

                    sendMessageToSocket(winPosString, identityString, winPositionSocket);

                    simsReady = false;

                    simsRunning = true;
                }
            }
            else
            {
                std::cout << "No simulations are waiting for initial window position" << std::endl;
            }

        }
        // send 20 jitter steps
        else if (input == 3)
        {
            if (simsRunning)
            {
                std::cout << "Send 20 jitter steps to the simulations" << std::endl;
                
                for (uint i = 0; i < 20; i++)
                {
                
                    std::cout << "Step counter: " << stepCounter << std::endl;
                
                    std::string jitterString =  std::to_string(stepVec.at(stepCounter)) + " "
                                                + std::to_string(yawVec.at(stepCounter)) + " "
                                                + std::to_string(pitchVec.at(stepCounter)) + " "
                                                + std::to_string(rollVec.at(stepCounter));
                
                    stepCounter++;
                
                    // send position string to identity on socket 3
                
                    sendMessageToSocket(jitterString, "", jitterSocket);
                }

            }
            else
            {
                std::cout << "No simulations are running" << std::endl;
            }

        }
        // change the window size
        else if (input == 4)
        {
            if (simsRunning)
            {
                for (uint i = 0; i < simNum; i++)
                {
        
                    std::cout << "Send the window position to identity: " << identityVec.at(i) << std::endl;
                            
                    std::string identityString = identityVec.at(i);
        
                    std::tuple<int, int, int, int, int, int, int, int, int, int> winPos = starPositionVec.at(0);
        
                    std::string winPosString =  std::to_string(std::get<5>(winPos)) + " " +
                                                std::to_string(std::get<6>(winPos)) + " " +
                                                std::to_string(std::get<7>(winPos)) + " " +
                                                std::to_string(std::get<8>(winPos)) + " " +
                                                std::to_string(std::get<9>(winPos));
        
                    std::cout << "Identity: " << identityString << std::endl;
        
                    std::cout << "Message: " << winPosString << std::endl;
        
                    // send position string to identity
        
                    sendMessageToSocket(winPosString, identityString, winPositionSocket);
        
                }
            }
            else
            {
                std::cout << "No simulations are running" << std::endl;
            }
        }
        // end running simulations
        else if (input == 5)
        {
            if (simsRunning || simsStarted)
            {
                simsStarted = false;

                std::cout << "End the simulations" << std::endl;

                for (uint i = 0; i < identityVec.size(); i++)
                {
                    sendMessageToSocket("6 6 0 0 0", identityVec.at(i), winPositionSocket);

                    sendMessageToSocket("", identityVec.at(i), jitterSocket);
                }

                simsRunning = false;
            }
            else
            {
                std::cout << "No simulations are running" << std::endl;
                
                continue;   
            }
        }
        //remove all hdf5 and txt files in the directory
        else if (input == 6)
        {
            std::string removeString = "rm ./*.hdf5 ./*.txt";

            const char* removeCall = removeString.c_str(); 

            system(removeCall);
        }
        else
        {
            std::cout << "Please enter a valid number between 0 and 6" << std::endl;
        }
    }
}











int main()
{

    // get data from jitter file

    //create an ifstream
    std::ifstream jitterFile(jitterFileName);

    double step, yaw, pitch, roll;

    // read the file line by line and save the values in the respective vector
    while (jitterFile >> step >> yaw >> pitch >> roll)
    {

        stepVec.push_back(step);
        yawVec.push_back(yaw);
        pitchVec.push_back(pitch);
        rollVec.push_back(roll);        
    }

    std::cout << "jitterVec size: " << stepVec.size() << std::endl;

    // get data from star position file

    int starId, subFieldRows1, subFieldCols1, subFieldRows2, subFieldCols2, rowPos1, rowPos2, colPos1, colPos2, orientation;

    //create an ifstream
    std::ifstream windowFile(winPositionFileName);

    // read all information from the file to a vector
    while (windowFile >> starId >> subFieldCols1 >> subFieldRows1 >> colPos1 >> rowPos1 >> subFieldCols2 >> subFieldRows2 >> colPos2 >> rowPos2 >> orientation)
    {
        starPositionVec.push_back(std::make_tuple(starId, subFieldRows1, subFieldCols1, rowPos1, colPos1, subFieldRows2, subFieldCols2, rowPos2, colPos2, orientation));
    }

    std::cout << "starPositionVec size: " << starPositionVec.size() << std::endl;


    // create the socket connections


    // define the context

    zmq::context_t context(1);

    // connect to platosim imagette socket (PUB-SUB pattern as subscriber)

    zmq::socket_t platoSimImagetteSocket (context, ZMQ_ROUTER);

    platoSimImagetteSocket.bind(platoSimImagetteAddress);

    // connect to platosim input socket (DEALER-CLIENT pattern as DEALER)

    zmq::socket_t platoSimWinPositionSocket (context, ZMQ_ROUTER);

    platoSimWinPositionSocket.bind(platoSimWinPositionAddress);

    // connect to platosim jitter socket (PUB-SUB pattern as publisher)

    zmq::socket_t platoSimJitterSocket (context, ZMQ_PUB);

    platoSimJitterSocket.bind(platoSimJitterAddress);    

    

    // start the sending and receiving thread

    std::thread receivingThread(receiveMessages, &platoSimImagetteSocket, &platoSimWinPositionSocket);

    std::thread sendingThread(sendMessages, &platoSimJitterSocket, &platoSimWinPositionSocket);


    // join the threads

    receivingThread.join();

    sendingThread.join();


    return 0;

}
