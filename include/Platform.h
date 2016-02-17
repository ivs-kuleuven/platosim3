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
#include "Telescope.h"
#include "ConfigurationParameters.h"
#include "JitterGenerator.h"


using namespace std;


class Telescope;  // forward declaration



class Platform : public Heartbeat, HDF5Writer
{
    public:

        Platform(ConfigurationParameters configParams, HDF5File &hdf5File, JitterGenerator &jitterGenerator);
        ~Platform();

        virtual void configure(ConfigurationParameters &configParams);

        void updatePointingCoordinates(Telescope const &telescope, double time);
        void setPointingCoordinates(double rightAscencsion, double declination, Unit unit = Angle::degrees);
        pair<double, double> getCurrentPointingCoordinates();

        virtual double getHeartbeatInterval() override;


    protected:

        arma::colvec rotateYawPitchRoll(arma::colvec coord, const double yaw, const double pitch, const double roll);
        

        double internalTime;                        // [s]
        double currentRA;                           // Current right Ascension of spacecraft pointing axis (zSC-axis) [rad]
        double currentDec;                          // Current declination     of spacecraft pointing axis (zSC-axis) [rad]
        double originalRA;                          // Original user-given right Ascension of spacecraft pointing axis (zSC-axis) [rad]
        double originalDec;                         // Original user-given declination     of spacecraft pointing axis (zSC-axis) [rad]

        JitterGenerator &jitterGenerator; 
 
    private:

};



#endif
