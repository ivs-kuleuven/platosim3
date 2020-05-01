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

    useStellarSpectra                  = configParams.getBoolean("SpectralDependency/GenerateFluxBasedOnTemperature");
    referenceWavelength                = configParams.getDouble("Camera/ThroughputLambdaC");
    binnumber                          = configParams.getInteger("SpectralDependency/NumberOfSpectralBins");
    lowerWavelength                    = configParams.getDouble("SpectralDependency/LowerWavelengthLimit");
    binwidth                           = configParams.getDouble("SpectralDependency/SpectralBinWidth");
    transmissionEfficiencySpectralBOL  = configParams.getDoubleVector("SpectralDependency/TransmissionEfficiencySpectralBOL");
    transmissionEfficiencySpectralEOL  = configParams.getDoubleVector("SpectralDependency/TransmissionEfficiencySpectralEOL");
    QESpectral                         = configParams.getDoubleVector("SpectralDependency/QuantumEfficiencySpectral");

    useQE                              = configParams.getBoolean("CCD/IncludeQuantumEfficiency");
    fileHasTemp                        = configParams.getBoolean("SpectralDependency/StarCatalogueHasTemp");

    meanQE                             = configParams.getDouble("CCD/QuantumEfficiency/MeanQuantumEfficiency");
    missionDuration                    = configParams.getDouble("ObservingParameters/MissionDuration") * 31536000.0; // [s]

    if ((transmissionEfficiencySpectralBOL.size() != binnumber) || (transmissionEfficiencySpectralEOL.size() != binnumber) || (QESpectral.size() != binnumber))
    {
        string errorMessage = "Spectral: number of transmission/QE values different from number of wavelength bins.";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
    }

    if (!fileHasTemp)
    {
        string errorMessage = "Spectral: StarCatalogue does not contain stellar temperature.";
        Log.warning(errorMessage);
    }

}



/**
 * \brief Return the transmission efficiency of the telescope in all spectral bins (number between 0 and 1)
 * 
 */

vector<double> SpectralDependenceUtility::getSpectralTransmissionEfficiency(double time)
{
    vector<double> transmission = transmissionEfficiencySpectralBOL;
    for (int i=0; i<transmissionEfficiencySpectralBOL.size(); i++)
    {
       transmission[i]-= (transmissionEfficiencySpectralBOL[i] - transmissionEfficiencySpectralEOL[i]) / missionDuration * time;
    }
    return transmission;
}

