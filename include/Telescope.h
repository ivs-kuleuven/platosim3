#ifndef TELESCOPE_H
#define TELESCOPE_H

#include <string>
#include <algorithm>

#include "Logger.h"
#include "Units.h"
#include "Constants.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "Platform.h"
#include "DriftGenerator.h"

using namespace std;


class Platform;  // Forward declaration


class Telescope  : public Heartbeat, HDF5Writer
{
    
    public:

        Telescope(ConfigurationParameters &configParams, HDF5File &hdf5File, Platform &platform, DriftGenerator &driftGenerator);
        virtual ~Telescope();

        virtual void configure(ConfigurationParameters &configParam);

        virtual void updateParameters(double time);

        arma::mat getPlatformToDriftedTelescopeRotationMatrix();
        arma::mat getDriftedTelescopeToPlatformRotationMatrix();

        arma::mat getPlatformToUndriftedTelescopeRotationMatrix();
        arma::mat getUndriftedTelescopeToPlatformRotationMatrix();


        virtual double getHeartbeatInterval() override;
        virtual double getTransmissionEfficiency(double time);
        
        double getLightCollectingArea();


    protected:

        virtual void initHDF5Groups() override;
        virtual void flushOutput() override;

        double originalAzimuthAngle;           // Original azimuth angle of telescope on platform in the inputfile  [rad]
        double originalTiltAngle;              // Original tilt angle of telescope on platform in the inputfile     [rad]
        double currentAlphaOpticalAxis;        // Current right ascension of the optical axis                       [rad]
        double currentDeltaOpticalAxis;        // Current declination of the optical axis                           [rad]
        double lightCollectingArea;            // Effective light collective area                                   [cm^2]
        double transmissionEfficiencyBOL;      // Efficiency at Beginning Of Life in [0,1]
        double transmissionEfficiencyEOL;      // Efficiency at End Of Life in [0,1]
        double missionDuration;                // Duration of the PLATO Mission, used for degrading parameters      [s]
        double driftYawRms;                    // RMS of thermo-elastic drift in yaw                                [rad]
        double driftPitchRms;                  // RMS of thermo-elastic drift in pitch                              [rad]
        double driftRollRms;                   // RMS of thermo-elastic drift in roll                               [rad]
        double driftTimeScale;                 // Timescale of thermo-elastic drift                                 [s]
 
        bool useDrift;                         // If false, the yaw, pitch, and roll of the thermo-elastic drift are always zero.
        bool writeTelescopeACS;                // Determines if the Transmission Efficiency should be writen to the HDF5 file

        double originalFocalPlaneOrientation;  // As in the input file [rad]

        vector<double> historyTime;            // The following vectors stores the telescope orientation angles and pointings
        vector<double> historyRA;              //     to write in the HDF5 file
        vector<double> historyDec;
        vector<double> historyYaw; 
        vector<double> historyPitch;
        vector<double> historyRoll;

        arma::mat rotDriftedTelescopeToSpacecraft;  // rotation matrix 
        arma::mat rotSpacecraftToDriftedTelescope;  // rotation matrix 

        arma::mat getUndriftedToDriftedRotationMatrix(const double yaw, const double pitch, const double roll);
        virtual void updateTelescopeOrientation(double time);

    private:

        double internalTime;               // Internal clock

        DriftGenerator &driftGenerator; 
        Platform &platform;
};

#endif
