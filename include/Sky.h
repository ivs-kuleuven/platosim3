
#ifndef SKY_H
#define SKY_H

#include <string>
#include <vector>

#include "Logger.h"
#include "TimeTicker.h"
#include "HDF5Writer.h"
#include "ConfigurationParameters.h"
#include "StarCatalog.h"

using namespace std;




class Sky : public TimeTicker, Hdf5Writer
{
    public:

        Sky(ConfigurationParameters configurationParameters);
        ~Sky();

        StarCatalog getStarsWithinRadiusFrom(double RA, double dec, double radius);
        double getSkyBackgroundFlux(double RA, double dec);
        double getZodiacalBackgroundFlux(double RA, double dec);
        double getUnresolvedStarsBackgroundFlux(double RA, double dec);

        // magnitude range
        

    protected:


    private:

};



#endif