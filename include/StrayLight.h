#ifndef STRAYLIGHT_H
#define STRAYLIGHT_H

#include <array>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "Camera.h"
#include "ConfigurationParameters.h"
#include "HDF5Writer.h"

// #include "HDF5File.h"
#include "armadillo"


class Camera;


struct GridPoint
{
    arma::vec point;
    double size;
};


class StrayLight : public HDF5Writer
{
  public:
    StrayLight(ConfigurationParameters &configParam, HDF5File &hdf5File,
               Camera &camera);
    // virtual ~StrayLight();

    void configure(ConfigurationParameters &configParam);
    // void updateParameters(double time);

  protected:
    std::vector<GridPoint> getGrid(double radius, unsigned int nPoints);
    void readInFile(std::string orbitPath, std::vector<arma::vec> &sc,
                    std::vector<arma::vec> &moon, std::vector<arma::vec> &sun);
    std::vector<std::string> splitLine(std::string &line);
    std::vector<arma::vec>
    getCelestialObjectGridSpectralRadiance(arma::vec sun, arma::vec object,
                                           double reflexivity,
                                           std::vector<GridPoint> &grid);
    std::array<double, 29> SolarSpectralIrradiance(double distance);
    std::tuple<std::vector<double>, std::vector<double>,
                std::vector<arma::vec>, std::vector<arma::vec>>
    getIrradianceAtCamera(Camera &camera, std::vector<GridPoint> grid,
                          std::vector<arma::vec> emmitterIrradiance,
                          arma::vec emmitterPosition, arma::vec cameraPosition);
    std::array<std::vector<std::array<double, 29>>, 5> getStrayLightAtDetector(
        std::array<std::vector<int>, 5> &rho_a,
        std::array<std::vector<std::array<double, 29>>, 5> &PST,
        std::vector<double> irradiance_alpha);
    std::array<std::vector<double>, 5>
    getNumberOfStraylightPhotoelectronsAtDetector(
        std::array<std::vector<std::array<double, 29>>, 5> &PST);
    std::array<double, 29> interpolatePST(double wl[3], double pst[3]);
    std::vector<double>
    extrapolate(std::vector<double> &irradiance_alpha, std::vector<int> &rho_AZ,
						       std::vector<std::array<double, 4>> &parameters);
    std::array<double,5> integrateOverGrid(std::array<std::vector<double>, 5> &strayLight);
    template<std::size_t N>
    double integrate(std::array<double, N> y, std::array<double, N> x);

    std::array<double, 4> getCubicParameters(double x_0, double x_1,
                                             double y_0, double y_1,
                                             double d0, double d1);
    std::pair<std::array<std::vector<int>, 5>,
              std::array<std::vector<std::array<double, 29>>, 5>>
    getPST(std::string pstPath);
    std::vector<arma::vec> moon;
    std::vector<arma::vec> sc;
    std::vector<arma::vec> sun;

    double radiusFOV;
    double pixelSize;
    double numExposure;
    double beginExposures;
    double cycleTime;

    int AZs[5];

    Camera &camera;
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
