#include "SpectralDependenceUtility.h"



/**
 * \brief constructor of the SpectralDependencyUtility class
 */
SpectralDependenceUtility::SpectralDependenceUtility(ConfigurationParameters &configParams)
{
    configure(configParams);
}



/**
 * \brief get needed variables from the input file
 *        
 */
void SpectralDependenceUtility::configure(ConfigurationParameters &configParams)
{

    useStellarSpectra               = configParams.getBoolean("SpectralDependency/GenerateFluxBasedOnTemperature");
    referenceWavelength             = configParams.getDouble("Camera/ThroughputLambdaC");
    binnumber                       = configParams.getInteger("SpectralDependency/NumberOfSpectralBins");
    lowerWavelength                 = configParams.getDouble("SpectralDependency/LowerWavelengthLimit");
    binwidth                        = configParams.getDouble("SpectralDependency/SpectralBinWidth");
    transmissionEfficiencySpectral  = configParams.getDoubleVector("SpectralDependency/TransmissionEfficiencySpectral");

    if (transmissionEfficiencySpectral.size() != binnumber)
    {
        string errorMessage = "Spectral: number of transmission values different from number of wavelength bins.";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
    }

}

