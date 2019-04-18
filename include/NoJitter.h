#ifndef NOJITTER_H
#define NOJITTER_H

#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>
#include <limits>

#include "Logger.h"
#include "Units.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "JitterGenerator.h"


using namespace std;



// Null Object Pattern for a null generator, always returning (yaw, pitch, roll)=(0,0,0)

class NoJitter: public JitterGenerator
{
    public:

    	static JitterGenerator* Instance();
        
        ~NoJitter();

        virtual void configure(ConfigurationParameters &configParams);
        virtual tuple<double, double, double> getNextYawPitchRoll(double time) override;
        virtual double getHeartbeatInterval() override;

    protected:

    	NoJitter();

    private:

 };



#endif
 
