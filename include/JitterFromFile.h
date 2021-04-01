#ifndef JITTERFROMFILE_H
#define JITTERFROMFILE_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

#include "Exceptions.h"
#include "Logger.h"
#include "Units.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "FileUtilities.h"
#include "JitterGenerator.h"


using namespace std;



class JitterFromFile : public JitterGenerator
{
    public:

        JitterFromFile(ConfigurationParameters &configurationParameters);
        ~JitterFromFile();

        virtual void configure(ConfigurationParameters &configParams);
        virtual tuple<double, double, double> getNextYawPitchRoll(double timeInterval) override;
        virtual double getHeartbeatInterval() override;

    protected:

        string pathToJitterFile;
        double beginTime;                 // Only read the jitter file from beginTime to endTime
        double endTime;

        int timeIndex;
        double internalTime;

        vector<double> timeFromFile;  // [s]
        vector<double> yaw;           // [rad]
        vector<double> pitch;         // [rad]
        vector<double> roll;          // [rad]


    private:

 };



#endif
 
