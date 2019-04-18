#ifndef NODRIFT_H
#define NODRIFT_H 

#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>

#include "Logger.h"
#include "Units.h"
#include "Exceptions.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "DriftGenerator.h"


using namespace std;



class NoDrift: public DriftGenerator
{
    public:

    	static DriftGenerator* Instance();

        ~NoDrift();

        virtual void configure();
        virtual tuple<double, double, double> getNextYawPitchRoll(double time) override;
        virtual double getHeartbeatInterval() override;

    protected:

        NoDrift();

    private:

 };



#endif
 
