#ifndef DRIFTGENERATOR_H
#define DRIFTGENERATOR_H

#include <string>
#include <vector>

#include "Logger.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"


using namespace std;



class DriftGenerator : public Heartbeat
{
    public:

        virtual ~DriftGenerator(){};

        virtual tuple<double, double, double> getNextYawPitchRoll(double time) = 0;

    protected:

    	DriftGenerator(){};

		static DriftGenerator* _instance;

    private:

 };



#endif