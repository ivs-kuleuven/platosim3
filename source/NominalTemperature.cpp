#include "NominalTemperature.h"



/**
 * \brief Constructor.
 *
 * \param configParams Configuration parameters from the input parameters file.
 *
 * \param component Component for which the temperature is fixed at the nominal operating temperature (FEE or CCD).
 */

NominalTemperature::NominalTemperature(ConfigurationParameters &configParams, string component)
{
	// Set the configuration parameters

	configure(configParams, component);
}










/**
 * \brief Destructor.
 */

NominalTemperature::~NominalTemperature()
{

}










/**
 * \brief Configure the temperature variations for the given component, based on the given configuration parameters.
 *
 * \param configParams Configuration parameters.
 *
 * \param component Component for which the temperature is fixed at the nominal operating temperature of the given component.
 */

void NominalTemperature::configure(ConfigurationParameters &configParams, string component)
{
    nominalTemperature = configParams.getDouble(component + "/NominalOperatingTemperature");
}










/**
 * \brief Get the nominal operating temperature.
 *
 * \param time  Time point for which the temperature is requested [s].
 *
 * \return Temperature [K].
 */
double NominalTemperature::getNextTemperature(double time)
{
	internalTime = time;

    return nominalTemperature;
}
