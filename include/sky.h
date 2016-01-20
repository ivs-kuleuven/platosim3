
#ifndef SKY_H
#define SKY_H

#include <string>
#include <vector>
#include "logger.h"
#include "timeticker.h"
#include "hdf5writer.h"
#include "configurationparameters.h"
#include "starcatalog.h"

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