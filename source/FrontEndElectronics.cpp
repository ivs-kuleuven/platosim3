#include "FrontEndElectronics.h"

FrontEndElectronics::FrontEndElectronics(ConfigurationParameters &configParam, HDF5File &hdf5file)
: HDF5Writer(hdf5file)
{
	// Parse the parameters from the configuration file.

	configure(configParam);

	// Generate FEE gain

	generateGain();
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
	nominalOperatingTemperature = configParam.getDouble("FEE/NominalOperatingTemp");

	readoutNoise  = configParam.getDouble("FEE/ReadoutNoise");

	refValueGain  = configParam.getDouble("FEE/Gain/RefValue");
	gainStability = configParam.getDouble("FEE/Gain/Stability");
	gainDelta     = configParam.getDouble("FEE/Gain/Delta");
	gainSeed      = configParam.getLong("RandomSeeds/FeeGainSeed");

	refValueBias  = configParam.getInteger("FEE/ElectronicOffset/RefValue");
	biasStability = configParam.getDouble("FEE/ElectronicOffset/Stability");
 }










/**
 * Generate FEE gain for the left and the right ADC.
 */

void FrontEndElectronics::generateGain()
{
	mt19937 gainGenerator;
	gainGenerator.seed(gainSeed);

	normal_distribution<double> gainDistribution = normal_distribution<double>(refValueGain, gainDelta * refValueGain);

	refValueGainLeft  = gainDistribution(gainGenerator);
	refValueGainRight = gainDistribution(gainGenerator);
}










/**
 * Returns the gain of the left ACD for the current operating temperature of the FEE.
 */

double FrontEndElectronics::getGainLeftAdc()
{
	return refValueGainLeft + gainStability * (getTemperature() - nominalOperatingTemperature);
}









/**
 * Returns the gain of the right ACD for the current operating temperature of the FEE.
 */

double FrontEndElectronics::getGainRightAdc()
{
	return refValueGainRight + gainStability * (getTemperature() - nominalOperatingTemperature);
}











/**
 * Returns the electronic offset for the current operating temperature of the FEE.
 */

double FrontEndElectronics::getElectronicOffset()
{
	return refValueBias + biasStability * (getTemperature() - nominalOperatingTemperature);
}











/**
 * \brief Returns the current temperature of the FEE, expressed in K.
 *
 * \note In the current implementation, we assume the temperature to be constant and equal to
 *       the nominal operating temperature of the FEE.
 */

double FrontEndElectronics::getTemperature()
{
	return nominalOperatingTemperature;
}
