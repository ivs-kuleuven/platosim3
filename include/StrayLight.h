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
class Detector;

struct GridPoint
{
    arma::vec point;
    double size;
};

struct CelestialObject
{
    double radius;
    double reflectivity;
};


class StrayLight : public HDF5Writer
{
  public:
    StrayLight(ConfigurationParameters &configParam, HDF5File &hdf5File,
               Camera &camera, Detector &detector);
    // virtual ~StrayLight();

    void configure(ConfigurationParameters &configParam);

    double getStrayLightMoon(double time);
    // void updateParameters(double time);
    double getStraylightFromAZ(const std::array<double, 5> &electronsAtDetector, double az);
protected:
    CelestialObject moon;
    CelestialObject earth;
  double getStrayLightObject(CelestialObject object, arma::vec sun_pos, arma::vec object_pos, arma::vec sc_pos, arma::Mat<double> telescopeAxis, unsigned int nGridPoints);


    std::vector<GridPoint> getGrid(double radius, unsigned int nPoints);
    void readInFile(std::string orbitPath, std::vector<arma::vec> &sc_pos,
                    std::vector<arma::vec> &moon_pos, std::vector<arma::vec> &sun_pos);
    std::vector<std::string> splitLine(std::string &line);
    std::vector<arma::vec>
    getCelestialObjectGridSpectralRadiance(arma::vec sun, arma::vec object,
                                           double reflectivity,
                                           std::vector<GridPoint> &grid);
    std::array<double, 29> solarSpectralIrradiance(double distance);
    std::tuple<std::vector<double>, std::vector<double>, std::vector<arma::vec>,
               std::vector<arma::vec>>
    getIrradianceAtCamera(Camera &camera, double row, double column,
                          std::vector<GridPoint> grid,
                          std::vector<arma::vec> emmitterIrradiance,
                          arma::vec emmitterPosition, arma::vec cameraPosition);
    std::array<std::vector<arma::vec>, 5> interpolatePSToverRho(
        std::array<std::vector<int>, 5> &rho_a,
        std::array<std::vector<std::array<double, 29>>, 5> &PST,
        std::vector<double> irradiance_alpha);
    std::array<std::vector<double>, 5>
    getNumberOfStraylightPhotoelectronsAtDetector(
        std::array<std::vector<arma::vec>, 5> &PST);
    std::array<double, 29> interpolatePST(double wl[3], double pst[3]);
    void extrapolatePST(std::array<std::vector<double>, 5> rho, std::array<std::vector<double>, 5> pst);
    double getPSTValue(double declination, double azimuth);
    std::vector<double> extrapolate(std::vector<double> &irradiance_alpha, std::vector<int> &rho, std::vector<std::array<double, 4>> &parameters);

    std::array<std::vector<std::array<double, 4>>,5> parameters;
    std::array<double, 4> getCubicParameters(double x_0, double x_1,
                                             double y_0, double y_1,
                                             double d0, double d1);

    void getPST(std::string pstPath);
    std::array<std::vector<double>, 5> rhoValues;
    std::vector<arma::vec> moon_positions;
    std::vector<arma::vec> sc_positions;
    std::vector<arma::vec> sun_positions;
    std::vector<double> times;


    std::array<int, 5> azs = {0, 45, 90, 135, 180};

    HDF5File pstFile;

    double radiusFOV;
    double pixelSize;
    double numExposure;
    double beginExposures;
    double cycleTime;
    arma::Mat<double> telescopeAxis{arma::Mat<double>(3,3)};

    int AZs[5];

    Camera &camera;
    Detector &detector;
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
