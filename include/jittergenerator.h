#ifndef JITTER_H
#define JITTER_H

#include <string>
#include <vector>

#include "Logger.h"
#include "ConfigurationParameters.h"


using namespace std;



class JitterGenerator
{
    public:

        JitterGenerator(ConfigurationParameters &configurationParameters);
        ~JitterGenerator();

        virtual void getNextYawPitchRoll(double &yaw, double &pitch, double &roll, double timeInterval);

    protected:


    private:

 };



#endif