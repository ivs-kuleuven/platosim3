#include "FrontEndElectronics.h"

FrontEndElectronics::FrontEndElectronics(ConfigurationParameters &configParam, HDF5File &hdf5file, TemperatureGenerator &temperatureGenerator)
: HDF5Writer(hdf5file), temperatureGenerator(temperatureGenerator)
{
	// Parse the parameters from the configuration file.

	configure(configParam);

	// Check gain

	checkGain();
}








/**
 * Destructor.
 *
 */

FrontEndElectronics::~FrontEndElectronics()
{
	flushOutput();
}










void FrontEndElectronics::configure(ConfigurationParameters &configParam)
 {
	nominalOperatingTemperature = configParam.getDouble("FEE/NominalOperatingTemperature");

	readoutNoise  = configParam.getDouble("FEE/ReadoutNoise");

	refValueGainLeft          = configParam.getDouble("FEE/Gain/RefValueLeft");
	refValueGainRight         = configParam.getDouble("FEE/Gain/RefValueRight");
	gainStability             = configParam.getDouble("FEE/Gain/Stability");
	gainAllowedDifference     = configParam.getDouble("FEE/Gain/AllowedDifference");

	refValueBias  = configParam.getInteger("FEE/ElectronicOffset/RefValue");
	biasStability = configParam.getDouble("FEE/ElectronicOffset/Stability");
 }







void FrontEndElectronics::checkGain()
{
	double allowedDifference = min(refValueGainLeft, refValueGainRight) * gainAllowedDifference / 100.0;

	if(abs(refValueGainLeft - refValueGainRight) > allowedDifference)
	{
		Log.warning("FrontEndElectornics: Difference in gain between ADC1 and ADC2 too large.");
	}
}









/**
 * Returns the gain of the left ACD for the current operating temperature of the FEE.
 */

double FrontEndElectronics::getGainLeftAdc(double time)
{
	return refValueGainLeft + gainStability * (getTemperature(time) - nominalOperatingTemperature);
}









/**
 * Returns the gain of the right ACD for the current operating temperature of the FEE.
 */

double FrontEndElectronics::getGainRightAdc(double time)
{
	return refValueGainRight + gainStability * (getTemperature(time) - nominalOperatingTemperature);
}











/**
 * Returns the electronic offset for the current operating temperature of the FEE.
 */

double FrontEndElectronics::getElectronicOffset(double time)
{
	return refValueBias + biasStability * (getTemperature(time) - nominalOperatingTemperature);
}











/**
 * \brief Returns the current temperature of the FEE, expressed in K.
 *
 * \note In the current implementation, we assume the temperature to be constant and equal to
 *       the nominal operating temperature of the FEE.
 */

double FrontEndElectronics::getTemperature(double time)
{
	return temperatureGenerator.getNextTemperature(time);
}
