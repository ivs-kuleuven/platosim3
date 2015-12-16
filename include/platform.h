#ifndef PLATFORM_H
#define PLATFORM_H

#include <string>
#include <vector>
#include "timeticker.h"
#include "hdf5writer.h"
#include "configurationparameters.h"
#include "jittergenerator.h"


using namespace std;




class Platform : public TimeTicker, Hdf5Writer
{
    public:

        Platform(ConfigurationParameters configurationParameters);
        ~Platform();

        void setPointingCoordinates(double rightAscencsion, double declination, double time);
        void updatePointingCoordinates(double &rightAscencsion, double &declination, double time);

    protected:

        double currentTime;                        // [s]
        double currentRA;                          // Right Ascension of pointing axis [rad]
        double currentDec;                         // Declination of pointing axis [rad]
        JitterGenerator* jitterGenerator; 
 
    private:

};



#endif
