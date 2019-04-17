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

	virtual void setCurrentJitterStep(bool endSimulation, double timeStep, double yaw, double pitch, double roll){};


    protected:


    private:

 };



#endif