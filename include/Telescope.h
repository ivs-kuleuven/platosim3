#ifndef TELESCOPE_H
#define TELESCOPE_H

#include <string>

#include "Logger.h"
#include "Units.h"
#include "Heartbeat.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"

using namespace std;




class Telescope  : public Heartbeat, HDF5Writer
{
	
	public:

		Telescope(ConfigurationParameters &configParams, HDF5File &hdf5File);
		~Telescope();

		virtual void configure(ConfigurationParameters &configParam);

		virtual void updatePointingCoordinates(double timeInterval);
		pair<double, double> getPointingCoordinates();

		double getTransmissionEfficiency();
		double getLightCollectingArea();
		double getFOVsolidAngle();

	protected:

		double alphaOpticalAxis;           // Current pointing right ascension [rad]
		double deltaOpticalAxis;           // Current pointing declination     [rad]
		double FOVradius;                  // Radius of the Field-of-view      [rad]
		double lightCollectingArea;        // Effective light collective area  [m^2]
		double transmissionEfficiency;     // in [0,1]
		double driftYawRms;                // RMS of thermo-elastic drift in yaw   [arcsec]
    	double driftPitchRms;              // RMS of thermo-elastic drift in pitch [arcsec]
    	double driftRollRms;               // RMS of thermo-elastic drift in roll  [arcsec]
    	double driftTimeScale;             // Timescale of thermo-elastic drift [s]

	private:

		//Platform platform;
};

#endif
