
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
#include "SkyCoordinates.h"
#include "TabulatedFunction.h"
#include "Skydata.h"
#include "StarCatalog.h"
#include "ConfigurationParameters.h"



using namespace std;




class Sky
{
    public:

        Sky(ConfigurationParameters &configParams);
        virtual ~Sky();

        virtual void configure(ConfigurationParameters &configParams);

        StarCatalog getStarsWithinRadiusFrom(double RA, double dec, double radius, Unit angleUnit = Angle::degrees);

        double zodiacalFlux(double RA, double dec, double lambda1, double lambda2);
        double zodiacalFlux(double RA, double dec, vector<double> &lambda, vector<double> &throughput);
        double stellarBackgroundFlux(double RA, double dec, double lambda1, double lambda2);
        double stellarBackgroundFlux(double RA, double decl, vector<double> &lambda, vector<double> &throughput);

        pair<double, double> getCoordinatesOfStarWithID(int id, Unit outputAngleUnit = Angle::degrees);
        double getVmagnitudeOfStarWithID(int id);

    protected:

        double solarRadiantFlux(double lambda);
        double solarRadiantFlux(double lambda1, double lambda2);
        double solarRadiantFlux(vector<double> &lambda, vector<double> &throughput);

    private:

        StarCatalog starCatalog;
        string starInputfile;

        vector<double> integrand;
        TabulatedFunction<vector<double>> tabfunction;

        void locate(double x, const double *array, int N, int &index);

};



#endif
