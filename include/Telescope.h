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
		~Telescope();

		virtual void configure(ConfigurationParameters &configParam);

		virtual void updatePointingCoordinates(double time);
		pair<double, double> getCurrentPointingCoordinates();

		double getTransmissionEfficiency();
		double getLightCollectingArea();
		double getFOVsolidAngle();


		pair<double, double> platformToTelescopePointingCoordinates(double alphaPlatfrom, double deltaPlatform);

		tuple<double, double, double> spacecraftToFocalPlaneCoordinates(const double xSC, const double ySC, const double zSC);
		tuple<double, double, double> focalPlaneToSpacecraftCoordinates(const double xFP, const double yFP, const double zFP);


	protected:

		double currentAlphaOpticalAxis;      // Current right ascension of the optical axis  [rad]
		double currentDeltaOpticalAxis;      // Current declination of the optical axis      [rad]
		double FOVsolidAngle;                // Solid angle of FOV of 1 telescope            [sr]
		double lightCollectingArea;          // Effective light collective area              [cm^2]
		double transmissionEfficiency;       // in [0,1]
		double driftYawRms;                  // RMS of thermo-elastic drift in yaw           [rad]
    	double driftPitchRms;                // RMS of thermo-elastic drift in pitch         [rad]
    	double driftRollRms;                 // RMS of thermo-elastic drift in roll          [rad]
    	double driftTimeScale;               // Timescale of thermo-elastic drift            [s]

	private:

		double internalTime;               // Internal clock
		Platform &platform;
};

#endif
