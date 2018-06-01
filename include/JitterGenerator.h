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

        virtual ~JitterGenerator(){};

        virtual tuple<double, double, double> getNextYawPitchRoll(double time) = 0;

    protected:
    	JitterGenerator(){};

      	static JitterGenerator* _instance;

    private:

 };




#endif