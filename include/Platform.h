#ifndef PLATFORM_H
#define PLATFORM_H

#include <string>
#include <vector>

#include "armadillo"

#include "Logger.h"
#include "Units.h"
#include "Constants.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "HDF5File.h"
#include "Telescope.h"
#include "ConfigurationParameters.h"
#include "JitterGenerator.h"


using namespace std;
using Constants::PI;


class Telescope;  // forward declaration



class Platform : public Heartbeat, HDF5Writer
{
    public:

        Platform(ConfigurationParameters configParams, HDF5File &hdf5File, JitterGenerator &jitterGenerator);
        virtual ~Platform();

        virtual void configure(ConfigurationParameters &configParams);

        void setPointingCoordinates(double rightAscencsion, double declination, Unit unit = Angle::degrees);
        pair<double, double> getPointingCoordinates(double time);
        pair<double, double> getInitialPointingCoordinates();


        virtual double getHeartbeatInterval() override;

        arma::colvec spacecraftToEquatorialCoordinates(arma::colvec &coordSC, bool useOriginalPointingCoordinates=false);




    protected:

        arma::colvec rotateYawPitchRoll(arma::colvec coord, const double yaw, const double pitch, const double roll);
        
        virtual void initHDF5Groups() override;
        virtual void flushOutput() override;


        bool useJitter;                             // If false, the yaw, pitch, and roll, are always zero.
        double internalTime;                        // [s]
        double currentRA;                           // Current right Ascension of spacecraft pointing axis (zSC-axis) [rad]
        double currentDec;                          // Current declination     of spacecraft pointing axis (zSC-axis) [rad]
        double originalRA;                          // Original user-given right Ascension of spacecraft pointing axis (zSC-axis) [rad]
        double originalDec;                         // Original user-given declination     of spacecraft pointing axis (zSC-axis) [rad]

        JitterGenerator &jitterGenerator; 
 
        vector<double> historyTime;                 // The following vectors stores all computed ASC values to write to HDF5
        vector<double> historyRA;
        vector<double> historyDec;
        vector<double> historyYaw;
        vector<double> historyPitch;
        vector<double> historyRoll;

    private:

};



#endif
