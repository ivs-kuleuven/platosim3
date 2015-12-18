#ifndef TELESCOPE_H
#define TELESCOPE_H

#include <string>
#include "timeticker.h"
#include "hdf5writer.h"
#include "platform.h"
#include "configurationparameters.h"

using namespace std;

class Telescope  : public TimeTicker, Hdf5Writer
{
	
	public:

		Telescope(ConfigurationParameters configurationParameters, Platform platform);
		~Telescope();

		Platform getPlatform();
		double getLightCollectingArea();
		double updatePointingCoordinates(double &alphaOpticalAxis, double &deltaOpticalAxis, double currentTime);

	protected:

		double alphaOpticalAxis;           // Current pointing right ascension [rad]
		double deltaOpticalAxis;           // Current pointing declination     [rad]

	private:

		Platform platform;
		void setPlatform(Platform platform);

		double lightCollectingArea;
		void setLightCollectingArea(double lightCollectingArea);

};

#endif
