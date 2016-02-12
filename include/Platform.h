#ifndef PLATFORM_H
#define PLATFORM_H

#include <string>
#include <vector>

#include "armadillo"

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

        Platform(ConfigurationParameters configParams, HDF5File &hdf5File, JitterGenerator &jitterGenerator);
        ~Platform();

        virtual void configure(ConfigurationParameters &configParams);

        void updatePointingCoordinates(double time);
        void setPointingCoordinates(double rightAscencsion, double declination, Unit unit = Angle::degrees);
        pair<double, double> getPointingCoordinates();

        virtual double getHeartbeatInterval() override;


    protected:

        arma::colvec rotateYawPitchRoll(arma::colvec coord, const double yaw, const double pitch, const double roll);
        

        double internalTime;                        // [s]
        double currentRA;                           // Right Ascension of pointing axis [rad]
        double currentDec;                          // Declination of pointing axis     [rad]

        JitterGenerator &jitterGenerator; 
 
    private:

};



#endif
