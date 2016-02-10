#ifndef JITTERGENERATOR_H
#define JITTERGENERATOR_H

#include <string>
#include <vector>

#include "Logger.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"


using namespace std;



class JitterGenerator : public Heartbeat
{
    public:

        JitterGenerator(){};
        ~JitterGenerator(){};

        virtual void getNextYawPitchRoll(double &yaw, double &pitch, double &roll, double timeInterval) = 0;

    protected:


    private:

 };



#endif