#ifndef FRONTENDELECTRONICS_H
#define FRONTENDELECTRONICS_H

#include <string>
#include <cmath>
#include <random>
#include <functional>

#include "armadillo"

#include "Constants.h"
#include "ArrayOperations.h"
#include "TemperatureGenerator.h"
#include "ConfigurationParameters.h"
#include "HDF5File.h"
#include "HDF5Writer.h"
#include "Logger.h"
#include "Units.h"


using namespace std;


class FrontEndElectronics: public HDF5Writer 
{
	public:

		FrontEndElectronics(ConfigurationParameters &configParam, HDF5File &hdf5File, TemperatureGenerator & temperatureGenerator);
		virtual ~FrontEndElectronics();

		virtual void configure(ConfigurationParameters &configParam);
		virtual void checkGain();

		virtual double getNominalOperatingTemperature(){return nominalOperatingTemperature;}	// Nominal operating temperature of the FEE [K]
		virtual double getTemperature(double time);                                			// FEE temperature at the given time

		virtual double getReadoutNoise(){return readoutNoise;};					                // FEE readout noise [e-/pixel]

		virtual double getGainLeftAdc(double time);
		virtual double getGainRightAdc(double time);

		virtual double getElectronicOffset(double time);

	protected:

		double nominalOperatingTemperature;		// Nominal operating temperature of the FEE [K]

		double readoutNoise;                     // FEE readout noise [e-/pixel]

		double refValueGainLeft;				 	// Reference value for the gain on the ACD reading the left-hand side of the detector [µV/e-]
		double refValueGainRight;				// Reference value for the gain on the ACD reading the right-hand side of the detector [µV/e-]
		double gainStability;					// Gain stability for the FEE [µV/e-/K]
		double gainAllowedDifference;			// Allowed difference in gain between both ADCs [% of the gain values]

		unsigned int refValueBias;				 // Reference value for the electronic offset [ADU/pixel]
		double biasStability;					 // Bias stability [ADU/pixel/K]

		TemperatureGenerator &temperatureGenerator;

	private:

};

#endif
