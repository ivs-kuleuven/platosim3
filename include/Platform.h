#ifndef PLATFORM_H
#define PLATFORM_H

#include <string>
#include <vector>

#include "Logger.h"
#include "Units.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "HDF5File.h"
#include "ConfigurationParameters.h"
#include "JitterGenerator.h"


using namespace std;




class Platform : public Heartbeat, HDF5Writer
{
    public:

        Platform(ConfigurationParameters configurationParameters, HDF5File &hdf5File, JitterGenerator &jitterGenerator);
        ~Platform();

        void updatePointingCoordinates(double time);
        void setPointingCoordinates(double rightAscencsion, double declination, Unit unit = Angle::degrees);
        pair<double, double> getPointingCoordinates();



    protected:

        double internalTime;                        // [s]
        
        double currentRA;                           // Right Ascension of pointing axis [rad]
        double currentDec;                          // Declination of pointing axis [rad]
        
        JitterGenerator &jitterGenerator; 
 
    private:

};



#endif
