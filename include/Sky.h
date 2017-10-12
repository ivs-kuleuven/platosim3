
#ifndef SKY_H
#define SKY_H

#include <string>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <valarray>
#include <tuple>

#include "Logger.h"
#include "Exceptions.h"
#include "Units.h"
#include "Constants.h"
#include "SkyCoordinates.h"
#include "TabulatedFunction.h"
#include "Skydata.h"
#include "Platform.h"
#include "ConfigurationParameters.h"



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
        tuple<unsigned int, double, double, double> getSelectedStar(unsigned int n);

        tuple<double, double, double> getInfoOfStarWithID(unsigned int starID);


        double zodiacalFlux(double RA, double dec, double lambda1, double lambda2);
        double zodiacalFlux(double RA, double dec, vector<double> &lambda, vector<double> &throughput);
        double stellarBackgroundFlux(double RA, double dec, double lambda1, double lambda2);
        double stellarBackgroundFlux(double RA, double decl, vector<double> &lambda, vector<double> &throughput);

        pair<double, double> getSunCoordinates(double julianDate, Unit outputAngleUnit = Angle::degrees);

    protected:

        double solarRadiantFlux(double lambda);
        double solarRadiantFlux(double lambda1, double lambda2);
        double solarRadiantFlux(vector<double> &lambda, vector<double> &throughput);

    private:

        map<unsigned int, tuple<double, double, double>> starDB;    // star database: starDB[stardID] contains (RA, dec, Vmag), 
                                                                    // with (RA, dec) in radians, and Vmag = Johnson V magnitude,
                                                                    // and starID the star identification number.
       
        vector<unsigned int> selectedStarID;  // Indices of the stars selected with selectStarsWithinRadiusFrom()
                                              // E.g. starID[selectedStarIndex[0]] is the ID of the first selected star 
        vector<double> selectedRA;            //      Corresponding selected Right Ascension [rad]
        vector<double> selectedDec;           //      Corresponding Declination              [rad]
        vector<double> selectedVmag;          //      Corresponding Johnson V magnitude

        string starInputfile;

        vector<double> integrand;
        TabulatedFunction<vector<double>> tabfunction;

        void locate(double x, const double *array, int N, int &index);

};



#endif
