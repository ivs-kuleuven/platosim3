#include"StrayLight.h"
#include "ConfigurationParameters.h"
#include "Constants.h"








double operator-(Time t1, Time t2);
double operator-(Time t1, Time t2)
{
    return difftime(t1.t, t2.t);
}



















StrayLight::StrayLight(ConfigurationParameters &configParam,
                     HDF5File &hdf5File,
		     Camera &camera)
: HDF5Writer(hdf5File), camera(camera)
{
    // Parse the parameters from the configuration file.

    configure(configParam);
}






void StrayLight::configure(ConfigurationParameters &configParam)
{
    numExposure		= configParam.getInteger("ObservingParameters/NumExposures");
    beginExposures	= configParam.getInteger("ObservingParameters/BeginExposureNr");
    cycleTime		= configParam.getInteger("ObservingParameters/CycleTime");
    radiusFOV		= deg2rad(configParam.getDouble("CCD/RelativeTransmissivity/RadiusFOV")); // [deg]

    double radiusMoon = 1.7381e6; // [m]
    double moon_reflectivity = 0.11;

    // Read in the positions of the file and save them into the vectors.

    std::string orbitPath =
        configParam.getAbsoluteFilename("StrayLight/FilePath");

    readInFile(orbitPath, sc, moon, sun);

    // Let's do the moon

    std::vector<GridPoint> grid;
    grid = getGrid(radiusMoon, 100);

    std::vector<arma::vec> celestialObjectSpectralRadiance =
        getCelestialObjectGridSpectralRadiance(sun[0], moon[0],
		moon_reflectivity, grid);

    getIrradianceAtCamera(camera, grid, celestialObjectSpectralRadiance, moon[0], sc[0]);
}
















/**
 * \brief: Calculates the spectral radiance for all gridpoints.
 *
 * \param: sun           Position of the sun.
 * \param: object        Position of the celestial object.
 * \param: reflexivity   Reflexivity of the celestial object.
 * \param: grid          The grid around the celstial object.
 *
 */
std::vector<arma::vec>
StrayLight::getCelestialObjectGridSpectralRadiance(arma::vec sun, arma::vec object,
                                       double reflexivity, std::vector<GridPoint> &grid)
{

    // Normalized vector from object to sun.
    arma::vec sunObejectVector = (sun - object);
    double distanceSunObject = sqrt(dot(sunObejectVector, sunObejectVector));
    sunObejectVector = sunObejectVector / distanceSunObject;

    // Solar irradiance at the object for different wavelengths
    std::array<double, 29> solarIrradiance = SolarSpectralIrradiance(distanceSunObject);
    std::vector<arma::vec> objectIrradiance;
    std::vector<GridPoint> new_grid;


    for (const GridPoint &p : grid)
    {
        // Normalize the vector from object to gridpoint.
        arma::vec x = p.point;
        x = x / sqrt(arma::as_scalar( x.t()*x));


        // Angle between sun, object and gridpoint.

        double cos_gamma = arma::as_scalar(x.t() * sunObejectVector);

        if (cos_gamma < 0)
        {
            // Dark side of the object
            continue;
        } else {
            // Bright side of the object

            // reflexivity of the celestial object
            double a = reflexivity * cos_gamma * cos_gamma / Constants::PI;

            // cosider the solar spectral irradiance using Plank's law
            arma::vec irradiance(29);
            for (int i = 0; i < 29; i++)
            {
                irradiance[i] = a * solarIrradiance[i];
            }
            objectIrradiance.push_back(irradiance);
            new_grid.push_back(p);
        }
    }
    grid = new_grid;
    return objectIrradiance;
}




















/**
 * \brief: This function returns a grid with nPoints points, around a celestial
 *         object with radius radius.
 *
 * \param: radius   Radius of the celstial object.
 * \param: nPoints  Amount of points in the grid.
 *
 * \output: returns a vector of gridpoints. (consisting of the position of the
 *          gridpoints and the surface area of the points.
 *
 */
std::vector<GridPoint> StrayLight::getGrid(double radius, int nPoints) {
    double dTheta = Constants::PI / nPoints;
    double dPhi   = 2*Constants::PI / nPoints;

    std::vector<GridPoint> grid;

    for (int i = 0; i < nPoints; i++)
    {
        for (int j = 0; j < nPoints; j++) {

            arma::vec gridPosition(3);
            gridPosition[0] = radius * sin(i * dTheta) * cos(j * dPhi);
            gridPosition[1] = radius * sin(i * dTheta) * sin( j * dPhi);
            gridPosition[2] = radius * cos(i * dTheta);

            GridPoint gridPoint;
            gridPoint.point = gridPosition;
            gridPoint.size = radius * radius * sin(i * dTheta) * dTheta * dPhi;

            grid.push_back(gridPoint);
        }
    }

    return grid;
}







std::vector<std::string> StrayLight::splitLine(std::string &line) {
    int i = 0;
    bool isString = false;
    std::vector<std::string> values;

    // Temporary string used to split string
    std::string s;
    while (line[i] != '\0') {
      if (isString) {
          if (line[i] == '\"') {
              isString = false;
          } else {
              s += line[i];
          }
      } else {
          if (line[i] == '\"') {
              isString = true;
          } else if (line[i] != ',') {
          // Append the char to the temp string
          s += line[i];
          } else {
          values.push_back(s);
          s.clear();
          }
      }
        i++;
    }
    return values;
}







void StrayLight::readInFile(std::string orbitPath, std::vector<arma::vec> &sc,
               std::vector<arma::vec> &moon, std::vector<arma::vec> &sun)
{
    Time t0 = Time("20260611T190026");

    double lower_bound = cycleTime * beginExposures;
    double upper_bound = cycleTime * numExposure + lower_bound;

    std::ifstream orbitFile(orbitPath);
    if (orbitFile.is_open())
    {
        std::string line;
        while (getline(orbitFile, line))
        {
            // Skip empty lines

            if (line.size() == 0) continue;

            // Skip lines that only contain white space

            const std::string whitespace = " /t/r/n";
            if (line.find_first_not_of(whitespace) == std::string::npos) continue;

            // Skip header line starting with '#'.

            if (line[0] == '#') continue;

            std::istringstream buffer(line);
            std::vector<std::string> value_of_line = splitLine(line);

            Time t = Time(value_of_line[0]);

            // Skip points that are not part of the simulation.

            if ( (t-t0 < lower_bound) || (t-t0 > upper_bound)) continue;

            arma::vec sc_row(3);
            sc_row[0] = std::stod(value_of_line[3]);
            sc_row[1] = std::stod(value_of_line[4]);
            sc_row[2] = std::stod(value_of_line[5]);

            arma::vec sun_row(3);
            sun_row[0] = std::stod(value_of_line[7]);
            sun_row[1] = std::stod(value_of_line[8]);
            sun_row[2] = std::stod(value_of_line[9]);

            arma::vec moon_row(3);
            moon_row[0] = std::stod(value_of_line[11]);
            moon_row[1] = std::stod(value_of_line[12]);
            moon_row[2] = std::stod(value_of_line[13]);

            sc.push_back(sc_row);
            sun.push_back(sun_row);
            moon.push_back(moon_row);
        }
    }
}







/**
 * \brief:    Calculates the solar irradiance using Planck's law
 *
 * \param:    distance  distance from object to the sun.
 *
 * \output:   Array containing the irradiance for different wavelengths
 *
 */
std::array<double, 29> StrayLight::SolarSpectralIrradiance(double distance) {

    // Wavelength range (400nm -> 1125nm)
    std::array<double, 29> irradiance;

    for (int i = 0; i < 29; i++)
    {
        double wavelength = 400 + i*25;
        double radiance = (2 * Constants::HPLANCK * Constants::CLIGHT * Constants::CLIGHT) /
                          (std::pow(wavelength, 5) *
                           (std::pow(Constants::E, (Constants::HPLANCK * Constants::CLIGHT / (wavelength * Constants::KBOLTZMANN * Constants::SOLART))) - 1));
        irradiance[i] = Constants::PI * radiance * std::pow((Constants::SOLARRADIUS / distance), 2);
    }
    return irradiance;
}






void StrayLight::getIrradianceAtCamera(Camera &camera, std::vector<GridPoint> grid, std::vector<arma::vec> emmitterIrradiance,  arma::vec emmitterPosition, arma::vec cameraPosition) {

    // Get the max angle where light can fall onto the camera

    double cos_alpha_max = cos(radiusFOV);


    // Camera pointing vectors

    double alpha, delta;
    double lambda, beta;

    tie(alpha, delta) = camera.focalPlaneToSkyCoordinates(0, 0);
    equatorial2ecliptic(alpha, delta, lambda, beta);

    arma::vec nCamera = {cos(lambda) * cos(beta), sin(lambda) * cos(beta), sin(beta) };


    for (int idx = 0; idx < grid.size(); idx++) {

        arma::vec n_Emmitter_Cam = emmitterPosition + grid[idx].point - cameraPosition;
        n_Emmitter_Cam = n_Emmitter_Cam / sqrt(arma::as_scalar(n_Emmitter_Cam.t() * n_Emmitter_Cam));

        double cos_alpha1 = arma::as_scalar(grid[idx].point.t() * n_Emmitter_Cam) / sqrt(arma::as_scalar(grid[idx].point.t() * grid[idx].point));
        double cos_alpha2 = arma::as_scalar(nCamera.t() * n_Emmitter_Cam);


        arma::vec E(29);

        if (cos_alpha1 < 0) {
            E.zeros();
        } else if (cos_alpha2 < 0) {
            E.zeros();
        } else {
            E = emmitterIrradiance[idx] * grid[idx].size * cos_alpha1;
        }

        double cos_gridAlpha = cos_alpha1;
        double cos_irraAlpha = cos_alpha2;

	arma::vec Eproj(29);
        if (cos_alpha2 < cos_alpha_max) {
            Eproj.zeros();
        }
        else {
	    Eproj = E*cos_alpha2;
        }
	// return E, Eproj, cos_irraAlpha, cos_gridAlpha
    }
    
}

Time::Time(std::string datetime)
{

    seconds = stoi(datetime.substr(13, 2));
    minutes = stoi(datetime.substr(11, 2));
    hours = stoi(datetime.substr(9, 2));
    days = stoi(datetime.substr(6, 2));
    months = stoi(datetime.substr(4, 2));
    years = stoi(datetime.substr(0, 4));

    struct tm tm;
    time_t rawtime;
    time(&rawtime);
    tm = *localtime(&rawtime);
    tm.tm_year = years - 1900;
    tm.tm_mon = months - 1;
    tm.tm_mday = days;
    tm.tm_hour = hours;
    tm.tm_min = minutes;
    tm.tm_sec = seconds;

    t = mktime(&tm);
}






    
