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

        virtual tuple<double, double, double> getNextYawPitchRoll(double time) = 0;

        virtual bool getSimulationState(){return endOfSimulation;};

    protected:

        bool endOfSimulation = false;

    private:

 };



#endif