#ifndef SPECTRALDEPENDENCEUTILITY_H
#define SPECTRALDEPENDENCEUTILITY_H

#include "ConfigurationParameters.h"
#include "armadillo"

// This class gets needed information to run PlatoSIm spectrally dependent

class SpectralDependenceUtility
{
    
    public:

        SpectralDependenceUtility(ConfigurationParameters &configParams);

        ~SpectralDependenceUtility() {};

        bool useStellarSpectra;
        double referenceWavelength;
        int binnumber;
        double lowerWavelength;
        double binwidth;
        vector<double> transmissionEfficiencySpectral;

    protected:


    private:

        virtual void configure(ConfigurationParameters &configParams);

};



#endif
