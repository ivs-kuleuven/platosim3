
#ifndef TELESCOPE_H
#define TELESCOPE_H

#include <string>
#include "platform.h"
#include "configurationparameters.h"


using namespace std;


class Telescope
{
    public:

        Telescope(ConfigurationParameters configurationParameters);
        ~Telescope(); 

 
    protected:

        double alphaOpticalAxis;                        // Current pointing right ascension [rad]
        double deltaOpticalAxis;                        // Current pointing declination [rad]
 
    private:
    
        Platform platform;


};



#endif
