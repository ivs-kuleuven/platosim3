#ifndef STRAYLIGHT_H
#define STRAYLIGHT_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <iostream>
#include <cmath>
#include <array>
#include <ctime>

#include "ConfigurationParameters.h"
#include "HDF5Writer.h"
// #include "HDF5File.h"
#include "armadillo"



struct GridPoint
{
    arma::vec point;
    double size;
};


class StrayLight : public HDF5Writer
{
    public:
        StrayLight(ConfigurationParameters &configParam, HDF5File &hdf5File);
        // virtual ~StrayLight();

        void configure(ConfigurationParameters &configParam);
       // void updateParameters(double time);

     protected:
         std::vector<GridPoint> getGrid(double radius, int nPoints);
         void readInFile(std::string orbitPath, std::vector<arma::vec> &sc,
                      std::vector<arma::vec> &moon,
                      std::vector<arma::vec> &sun);
         std::vector<std::string> splitLine(std::string &line);
         std::vector<std::array<double, 29>>
             getCelestialObjectGridSpectralRadiance(arma::vec sun, arma::vec object,
                                                     double reflexivity,
                                                     std::vector<GridPoint> &grid);
         std::array<double, 29> SolarSpectralIrradiance(double distance);

         std::vector<arma::vec> moon;
         std::vector<arma::vec> sc;
         std::vector<arma::vec> sun;

         double numExposure;
         double beginExposures;
         double cycleTime;
         };



class Time
{
    public:
        Time(std::string datetime);
            Time();

        time_t t;

private:
        int seconds;
        int minutes;
        int hours;
        int days;
        int months;
        int years;
};

#endif
