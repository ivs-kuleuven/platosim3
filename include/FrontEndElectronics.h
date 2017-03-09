#ifndef FRONTENDELECTRONICS_H
#define FRONTENDELECTRONICS_H

#include <string>
#include <cmath>
#include <random>
#include <functional>

#include "armadillo"

#include "Constants.h"
#include "ArrayOperations.h"
#include "ConfigurationParameters.h"
#include "HDF5File.h"
#include "HDF5Writer.h"
#include "Logger.h"
#include "Units.h"


using namespace std;


class FrontEndElectronics: public HDF5Writer 
{
	public:

		FrontEndElectronics(ConfigurationParameters &configParam, HDF5File &hdf5File);
		virtual ~FrontEndElectronics();

		virtual void configure(ConfigurationParameters &configParam);

		virtual double getNominalOperatingTemperature(){return nominalOperatingTemperature;}	// Nominal operating temperature of the FEE [K]
		virtual double getTemperature();						                                // Current FEE temperature

		virtual double getReadoutNoise(){return readoutNoise;};					                // FEE readout noise [e-/pixel]

		virtual double getGainLeftAdc();
		virtual double getGainRightAdc();

		virtual double getElectronicOffset();

	protected:

		virtual void generateGain();

		double nominalOperatingTemperature;		 // Nominal operating temperature of the FEE [K]

		double readoutNoise;                     // FEE readout noise [e-/pixel]

		double refValueGain;					 // Reference value for the FEE gain [µV/e-]
		double gainStability;					 // Gain stability for the FEE [µV/e-/K]
		double gainDelta;						 // Allowed difference in gain between the left and the right half of the detector [%]

		double refValueGainLeft;				 // Reference value for the gain on the ACD reading the left-hand side of the detector [µV/e-]
		double refValueGainRight;				 // Reference value for the gain on the ACD reading the right-hand side of the detector [µV/e-]

		unsigned int refValueBias;				 // Reference value for the electronic offset [ADU/pixel]
		double biasStability;					 // Bias stability [ADU/pixel/K]

		long gainSeed;

	private:

};

#endif
