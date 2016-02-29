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
        virtual ~JitterGenerator(){};

        virtual tuple<double, double, double> getNextYawPitchRoll(double timeInterval) = 0;

    protected:


    private:

 };



#endif