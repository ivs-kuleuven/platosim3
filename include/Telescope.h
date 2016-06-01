#ifndef TELESCOPE_H
#define TELESCOPE_H

#include <string>

#include "Logger.h"
#include "Units.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "Platform.h"

using namespace std;


class Platform;  // Forward declaration


class Telescope  : public Heartbeat, HDF5Writer
{
	
	public:

		Telescope(ConfigurationParameters &configParams, HDF5File &hdf5File, Platform &platform);
		virtual ~Telescope();

		virtual void configure(ConfigurationParameters &configParam);

		virtual void updatePointingCoordinates(double time);
        pair<double, double> getCurrentPointingCoordinates();
        pair<double, double> getInitialPointingCoordinates();

		double getTransmissionEfficiency();
		double getLightCollectingArea();


		pair<double, double> platformToTelescopePointingCoordinates(double alphaPlatfrom, double deltaPlatform);


	protected:

		virtual void initHDF5Groups() override;
		virtual void flushOutput() override;

		double azimuthAngle;                 // Azimuth angle of telescope on platform       [rad]
		double tiltAngle;                    // Tilt angle of telescope on platform          [rad]
		double currentAlphaOpticalAxis;      // Current right ascension of the optical axis  [rad]
		double currentDeltaOpticalAxis;      // Current declination of the optical axis      [rad]
		double lightCollectingArea;          // Effective light collective area              [cm^2]
		double transmissionEfficiency;       // in [0,1]
		double driftYawRms;                  // RMS of thermo-elastic drift in yaw           [rad]
    	double driftPitchRms;                // RMS of thermo-elastic drift in pitch         [rad]
    	double driftRollRms;                 // RMS of thermo-elastic drift in roll          [rad]
    	double driftTimeScale;               // Timescale of thermo-elastic drift            [s]


    	vector<double> historyTime;          // The following vectors stores the telescope pointings
        vector<double> historyRA;            //     to write in the HDF5 file
        vector<double> historyDec;

	private:

		double internalTime;               // Internal clock
		Platform &platform;
};

#endif
