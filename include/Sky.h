
#ifndef SKY_H
#define SKY_H

#include <string>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>

#include "Logger.h"
#include "Units.h"
#include "Constants.h"
#include "Coordinates.h"
#include "TabulatedFunction.h"
#include "Skydata.h"
#include "StarCatalog.h"
#include "ConfigurationParameters.h"



using namespace std;




class Sky
{
    public:

        Sky(ConfigurationParameters &configParams);
        ~Sky();

        virtual void configure(ConfigurationParameters &configParams);

        StarCatalog getStarsWithinRadiusFrom(double RA, double dec, double radius, Unit angleUnit = Angle::degrees);

        double ZodiacalFlux(const double RA, const double dec, const double lambda1, const double lambda2);
        double ZodiacalFlux(const double RA, const double dec, vector<double> &lambda, vector<double> &throughput);
        double StellarBackgroundFlux(const double RA, const double dec, const double lambda1, const double lambda2);
        double StellarBackgroundFlux(const double RA, const double decl, vector<double> &lambda, vector<double> &throughput);

    protected:

        double SolarRadiantFlux(const double lambda);
        double SolarRadiantFlux(const double lambda1, const double lambda2);
        double SolarRadiantFlux(vector<double> &lambda, vector<double> &throughput);

    private:

        StarCatalog starCatalog;
        string starInputfile;

        vector<double> integrand;
        TabulatedFunction<vector<double>> tabfunction;

        void locate(double x, const double *array, int N, int &index);

};



#endif
