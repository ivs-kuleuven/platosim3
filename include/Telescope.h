#ifndef TELESCOPE_H
#define TELESCOPE_H

#include <string>
#include <algorithm>

#include "Logger.h"
#include "Units.h"
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

		virtual void updatePointingCoordinates(double time);
        pair<double, double> getCurrentPointingCoordinates();
        pair<double, double> getInitialPointingCoordinates();

        virtual double getHeartbeatInterval() override;

		double getTransmissionEfficiency();
		double getLightCollectingArea();
        double getCurrentFocalPlaneOrientation();
        double getInitialFocalPlaneOrientation();


		pair<double, double> platformToTelescopePointingCoordinates(double alphaPlatfrom, double deltaPlatform);

		tuple<double, double, double> getNextYawPitchRoll(double timeInterval);
		tuple<double, double, double> getThermalDrift(double timeInterval);

	protected:

		virtual void initHDF5Groups() override;
		virtual void flushOutput() override;

		double originalAzimuthAngle;           // Original azimuth angle of telescope on platform in the inputfile  [rad]
		double originalTiltAngle;              // Original tilt angle of telescope on platform in the inputfile     [rad]
        double currentAzimuthAngle;            // Current azimuth angle of telescope on platform                    [rad]
        double currentTiltAngle;               // Current tilt angle of telescope on platform                       [rad]
		double currentAlphaOpticalAxis;        // Current right ascension of the optical axis                       [rad]
		double currentDeltaOpticalAxis;        // Current declination of the optical axis                           [rad]
		double lightCollectingArea;            // Effective light collective area                                   [cm^2]
		double transmissionEfficiency;         // in [0,1]
		double driftYawRms;                    // RMS of thermo-elastic drift in yaw                                [rad]
    	double driftPitchRms;                  // RMS of thermo-elastic drift in pitch                              [rad]
    	double driftRollRms;                   // RMS of thermo-elastic drift in roll                               [rad]
    	double driftTimeScale;                 // Timescale of thermo-elastic drift                                 [s]

        bool useDrift;                         // If false, the yaw, pitch, and roll of the thermo-elastic drift are always zero.

        double originalFocalPlaneOrientation;  // As in the input file [rad]
        double currentFocalPlaneOrientation;   // Perturbed due to thermo-elastic drift [rad]

    	vector<double> historyTime;            // The following vectors stores the telescope orientation angles and pointings
        vector<double> historyRA;              //     to write in the HDF5 file
        vector<double> historyDec;
        vector<double> historyYaw; 
        vector<double> historyPitch;
        vector<double> historyRoll;
        vector<double> historyAzimuth; 
        vector<double> historyTilt;
        vector<double> historyFocalPlaneOrientation;



	private:

		double internalTime;               // Internal clock

        DriftGenerator &driftGenerator; 
		Platform &platform;
};

#endif
