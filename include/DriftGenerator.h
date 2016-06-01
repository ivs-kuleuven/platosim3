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

        DriftGenerator(){};
        virtual ~DriftGenerator(){};

        virtual tuple<double, double, double> getNextYawPitchRoll(double timeInterval) = 0;

    protected:


    private:

 };



#endif