
#ifndef SKY_H
#define SKY_H

#include <string>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <valarray>
#include <memory>
#include <tuple>
#include <map>

#include "Logger.h"
#include "Exceptions.h"
#include "FileUtilities.h"
#include "StringUtilities.h"
#include "Units.h"
#include "Constants.h"
#include "SkyCoordinates.h"
#include "TabulatedFunction.h"
#include "Skydata.h"
#include "Platform.h"
#include "ConfigurationParameters.h"
#include "Parameter.h"



using namespace std;




class Sky
{
    public:

        Sky(ConfigurationParameters &configParams);
        virtual ~Sky();

        virtual void configure(ConfigurationParameters &configParams);
        virtual void updateParameters(double time);

        unsigned long selectStarsWithinRadiusFrom(double RA, double dec, double radius, Unit angleUnit = Angle::degrees);
        void aberrateSelectedStarPositions(Platform &platform, string aberrationCorrectionType, double startTime, double timeMiddle);
        tuple<unsigned int, double, double, double, double> getSelectedStar(unsigned int n);  //%% Added Temperature for spectral dependency

	unsigned long selectStarsWithinRadius();  //%% Added with spectral dependency, returns number of stars

        tuple<double, double, double, double> getInfoOfStarWithID(unsigned int starID);  //%% Added Temperature for spectral dependency


        double zodiacalFlux(double RA, double dec, double lambda1, double lambda2);
        double zodiacalFlux(double RA, double dec, vector<double> &lambda, vector<double> &throughput);
        double stellarBackgroundFlux(double RA, double dec, double lambda1, double lambda2);
        double stellarBackgroundFlux(double RA, double decl, vector<double> &lambda, vector<double> &throughput);

        pair<double, double> getSunCoordinates(double julianDate, Unit outputAngleUnit = Angle::degrees);

    protected:

        double solarRadiantFlux(double lambda);
        double solarRadiantFlux(double lambda1, double lambda2);
        double solarRadiantFlux(vector<double> &lambda, vector<double> &throughput);

        map<unsigned int, unique_ptr<Parameter<double>>> deltaMagnitude;     // deltaMagnitude[starID] contains pointer to time series file.  

    private:

        map<unsigned int, tuple<double, double, double, double>> starDB;    // star database: starDB[stardID] contains (RA, dec, Vmag), 
                                                                    // with (RA, dec) in radians, and Vmag = Johnson V magnitude,
                                                                    // and starID the star identification number.
                                                                     //%% Added Temperature for spectral dependency
       
        vector<unsigned int> selectedStarID;  // IDs of the stars selected with selectStarsWithinRadiusFrom()
        vector<double> selectedRA;            //      Corresponding selected Right Ascension [rad]
        vector<double> selectedDec;           //      Corresponding Declination              [rad]
        vector<double> selectedVmag;          //      Corresponding Johnson V magnitude
	vector<double> selectedTempStar;   //%% Added Temperature for spectral dependency

        vector<unsigned int> selectedVariableStars;    // Indices of those selected stars that have an entry in deltaMagnitude.

        string starInputfile;

        vector<double> integrand;
        TabulatedFunction<vector<double>> tabfunction;

        void locate(double x, const double *array, int N, int &index);

};



#endif
