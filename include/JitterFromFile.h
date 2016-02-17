#ifndef JITTERFROMFILE_H
#define JITTERFROMFILE_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

#include "Logger.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "JitterGenerator.h"


using namespace std;



class JitterFromFile : public JitterGenerator
{
    public:

        JitterFromFile(ConfigurationParameters &configurationParameters);
        ~JitterFromFile();

        virtual void configure(ConfigurationParameters &configParams);
        virtual void getNextYawPitchRoll(double &yaw, double &pitch, double &roll, double timeInterval) override;
        virtual double getHeartbeatInterval() override;

    protected:

        string pathToJitterFile;

        int timeIndex;

        vector<double> time;
        vector<double> yaw;      // [rad]
        vector<double> pitch;    // [rad]
        vector<double> roll;     // [rad]


    private:

 };



#endif
 