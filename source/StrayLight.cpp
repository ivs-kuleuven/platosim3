#include "StrayLight.h"
#include "ConfigurationParameters.h"
#include "Constants.h"
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

double operator-(Time t1, Time t2);
double operator-(Time t1, Time t2) { return difftime(t1.t, t2.t); }





/**
 * \brief: Constructor for StrayLight class.
 *
 * \param: configParam   Configuration parameters for the detector.
 * \param: hdf5file      HDF5 file to write to.
 * \param: camera        Camera on which the straylight will fall.
 *
 */
StrayLight::StrayLight(ConfigurationParameters &configParam, HDF5File &hdf5File,
                       Camera &camera)
    : HDF5Writer(hdf5File), camera(camera)
{
    // Parse the parameters from the configuration file.

    configure(configParam);
}





/**
 * \brief: Configure the StrayLight object using the configuration parameters.
 *
 * \param: configParam      Configuration parameters.
 *
 */
void StrayLight::configure(ConfigurationParameters &configParam)
{
    numExposure = configParam.getInteger("ObservingParameters/NumExposures");
    beginExposures =
        configParam.getInteger("ObservingParameters/BeginExposureNr");
    cycleTime = configParam.getInteger("ObservingParameters/CycleTime");
    radiusFOV = deg2rad(
        configParam.getDouble("CCD/RelativeTransmissivity/RadiusFOV")); // [deg]
    pixelSize = configParam.getDouble("CCD/PixelSize") * 1e-6;

    double radiusMoon = 1.7381e6; // [m]
    double moon_reflectivity = 0.11;

    // Read in the positions of the file and save them into the vectors.

    std::string orbitPath =
        configParam.getAbsoluteFilename("StrayLight/FilePath");
    readInFile(orbitPath, sc, moon, sun);

    // Read in the PST file and save it into vectors.

    std::string pstPath = configParam.getAbsoluteFilename("StrayLight/PstPath");
    std::array<std::vector<std::array<double, 29>>, 5> PST;
    std::array<std::vector<int>, 5> rho;

    tie(rho, PST) = getPST(pstPath);

    // Let's do the moon

    std::vector<GridPoint> grid;
    grid = getGrid(radiusMoon, 100);

    std::vector<arma::vec> celestialObjectSpectralRadiance =
        getCelestialObjectGridSpectralRadiance(sun[0], moon[0],
                                               moon_reflectivity, grid);


    std::vector<double> irradiance_alpha;
    std::vector<double> grid_alpha;
    std::vector<arma::vec> x;
    std::vector<arma::vec> y;
    std::tie(irradiance_alpha, grid_alpha, x, y) = getIrradianceAtCamera(
        camera, grid, celestialObjectSpectralRadiance, moon[0], sc[0]);


    std::array<std::vector<std::array<double, 29>>, 5> PST_interpolated = getStrayLightAtDetector(rho, PST, irradiance_alpha);
    std::array<std::vector<double>, 5> strayLightPhotoelectronsAtDetector =
        getNumberOfStraylightPhotoelectronsAtDetector(PST_interpolated);
    std::array<double,5> staylightAtDetector = integrateOverGrid(strayLightPhotoelectronsAtDetector);
}




std::array<double,5> StrayLight::integrateOverGrid(std::array<std::vector<double>, 5> &strayLight)
{
    std::cout << " " << std::endl;
    std::array<double, 5> totalElectrons;

    for (int i = 0; i < 5; i++)
    {
        totalElectrons[i] =
            std::accumulate(strayLight[i].begin(), strayLight[i].end(), 0);
    }
    std::cout << totalElectrons[0] << std::endl;
    return totalElectrons;
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
std::vector<arma::vec> StrayLight::getCelestialObjectGridSpectralRadiance(
    arma::vec sun, arma::vec object, double reflexivity,
    std::vector<GridPoint> &grid)
{

    // Normalized vector from object to sun.
    arma::vec sunObejectVector = (sun - object);
    double distanceSunObject = sqrt(dot(sunObejectVector, sunObejectVector));
    sunObejectVector = sunObejectVector / distanceSunObject;

    // Solar irradiance at the object for different wavelengths
    std::array<double, 29> solarIrradiance =
        SolarSpectralIrradiance(distanceSunObject);
    std::vector<arma::vec> objectIrradiance;
    std::vector<GridPoint> new_grid;

    for (const GridPoint &p : grid)
    {
        // Normalize the vector from object to gridpoint.
        arma::vec x = p.point;
        x = x / sqrt(arma::as_scalar(x.t() * x));

        // Angle between sun, object and gridpoint.

        double cos_gamma = arma::as_scalar(x.t() * sunObejectVector);

        if (cos_gamma < 0)
        {
            // Dark side of the object
            continue;
        }
        else
        {
            // Bright side of the object

            // reflexivity of the celestial object
            double a = reflexivity * cos_gamma * cos_gamma / Constants::PI;

            // cosider the solar spectral irradiance using Plank's law
            arma::vec irradiance(29);
            for (unsigned int i = 0; i < 29; i++)
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
std::vector<GridPoint> StrayLight::getGrid(double radius, unsigned int nPoints)
{
    double dTheta = Constants::PI / nPoints;
    double dPhi = 2 * Constants::PI / nPoints;

    std::vector<GridPoint> grid;

    // We should skip the case i=0, since this has size 0.

    for (unsigned int i = 1; i < nPoints; i++)
    {
        for (unsigned int j = 0; j < nPoints; j++)
        {

            arma::vec gridPosition(3);
            gridPosition[0] = radius * sin(i * dTheta) * cos(j * dPhi);
            gridPosition[1] = radius * sin(i * dTheta) * sin(j * dPhi);
            gridPosition[2] = radius * cos(i * dTheta);

            GridPoint gridPoint;
            gridPoint.point = gridPosition;
            gridPoint.size = radius * radius * sin(i * dTheta) * dTheta * dPhi;

            grid.push_back(gridPoint);
        }
    }

    return grid;
}





/**
 * \brief: Reads in a line from CSV file and returns a vector with elements the
 *         comma seperated variables.
 *
 * \param: line   The line from the CSV file. 
 */
std::vector<std::string> StrayLight::splitLine(std::string &line)
{
    int i = 0;
    bool isString = false;
    std::vector<std::string> values;

    // Temporary string used to split string

    std::string s;
    while (line[i] != '\0')
    {
        if (isString)
        {
            if (line[i] == '\"')
            {
                isString = false;
            }
            else
            {
                s += line[i];
            }
        }
        else
        {
            if (line[i] == '\"')
            {
                isString = true;
            }
            else if (line[i] != ',')
            {
                // Append the char to the temp string
                s += line[i];
            }
            else
            {
                values.push_back(s);
                s.clear();
            }
        }
        i++;
    }
    values.push_back(s);
    return values;
}






void StrayLight::readInFile(std::string orbitPath, std::vector<arma::vec> &sc,
                            std::vector<arma::vec> &moon,
                            std::vector<arma::vec> &sun)
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

            if (line.size() == 0)
                continue;

            // Skip lines that only contain white space

            const std::string whitespace = " /t/r/n";
            if (line.find_first_not_of(whitespace) == std::string::npos)
                continue;

            // Skip header line starting with '#'.

            if (line[0] == '#')
                continue;

            std::istringstream buffer(line);
            std::vector<std::string> value_of_line = splitLine(line);

            Time t = Time(value_of_line[0]);

            // Skip points that are not part of the simulation.

            if ((t - t0 < lower_bound) || (t - t0 > upper_bound))
                continue;

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

	    // We add these values and convert from [km] units to [m] units
            sc.push_back(sc_row * 1000);
            sun.push_back(sun_row * 1000);
            moon.push_back(moon_row * 1000);
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
std::array<double, 29> StrayLight::SolarSpectralIrradiance(double distance)
{

    // Wavelength range (400nm -> 1125nm)
    std::array<double, 29> irradiance;

    for (unsigned int i = 0; i < 29; i++)
    {
        double wavelength = 400 + i * 25;
        double radiance =
            (2 * Constants::HPLANCK * Constants::CLIGHT * Constants::CLIGHT) /
            (std::pow(wavelength, 5) *
             (std::pow(Constants::E, (Constants::HPLANCK * Constants::CLIGHT /
                                      (wavelength * Constants::KBOLTZMANN *
                                       Constants::SOLART))) -
              1));
        irradiance[i] = Constants::PI * radiance *
                        std::pow((Constants::SOLARRADIUS / distance), 2);
    }
    return irradiance;
}





/**
 * \brief Parses the pst file that is used as input and returns
 *        a matrix with interpolated pst values for a larger range
 *        of wavelengths.
 *
 * \details:
 *
 * Parses a PST file where we find for every AZ angle [0, 45, 90, 135, 180]
 * a matrix with pst values for the rho_a angle (in row) and the wavelength
 * (in column), but only for the wavelengths (500nm, 700nm and 1000nm).
 *
 * Return 5 PST matrices (for every AZ angle ) with every matrix has a row
 * for every angle rho_a and column for every wavelength (400nm -> 1125nm,
 * in steps of 25nm).
 * Also returns 5 vectors with the corresponding rho_a angle values.
 *
 * \param: pstPath  Path to the corresponding PST file.
 */
std::pair<std::array<std::vector<int>, 5>,
          std::array<std::vector<std::array<double, 29>>, 5>>
StrayLight::getPST(std::string pstPath)
{
    // Read PST file
    // Save PST file in datastructure

    std::array<std::vector<std::array<double, 4>>, 5> PSTs;
    std::ifstream pstFile(pstPath);

    if (pstFile.is_open())
    {

        int az_idx = 0;
        std::string line;
        std::vector<std::array<double, 4>> PST;
        while (getline(pstFile, line))
        {



            // Skip empty lines

            if (line.size() == 0)
                continue;

            // Skip lines that only contain white space

            const std::string whitespace = " /t/r/n";
            if (line.find_first_not_of(whitespace) == std::string::npos)
                continue;

            // Skip header line starting with '#'.

            if (line[0] == '#')
                continue;

            if (line.substr(0, 2) == "AZ")
            {
                int az = std::stoi(line.substr(3, 5));
                AZs[az_idx] = az;

                if (az != 0)
                {
                    PSTs[az_idx - 1] = PST;
                    PST.clear();
                }

                az_idx += 1;
            }
            else
            {
                std::istringstream buffer(line);
                std::vector<std::string> value_of_line = splitLine(line);

                std::array<double, 4> lineX;
                lineX[0] = std::stod(value_of_line[0]);
                lineX[1] = std::stod(value_of_line[1]);
                lineX[2] = std::stod(value_of_line[2]);
                lineX[3] = std::stod(value_of_line[3]);

                PST.push_back(lineX);
            }
        }
        PSTs[az_idx - 1] = PST;
    }

    std::array<std::vector<std::array<double, 29>>, 5> PST_interpolated;
    std::array<std::vector<int>, 5> rho_a;

    for (std::array<double, 4> pst : PSTs[1])
    {
        double pstx[3] = {pst[1], pst[2], pst[3]};
        double wl[3] = {500, 750, 1000};

        PST_interpolated[0].push_back(interpolatePST(wl, pstx));
        rho_a[0].push_back(int(pst[0]));
    }

    return std::make_pair(rho_a, PST_interpolated);
}





std::array<double, 29> StrayLight::interpolatePST(double wl[3], double pst[3])
{
    double N = (std::pow(wl[0], 2) - std::pow(wl[1], 2)) * (wl[0] - wl[2]) -
               (std::pow(wl[0], 2) - std::pow(wl[2], 2)) * (wl[0] - wl[1]);

    double a = (wl[0] - wl[2]) * (pst[0] - pst[1]) -
               (wl[0] - wl[1]) * (pst[0] - pst[2]);
    a = a / N;
    double b = (std::pow(wl[0], 2) - std::pow(wl[1], 2)) * (pst[0] - pst[2]) +
               (std::pow(wl[2], 2) - std::pow(wl[0], 2)) * (pst[0] - pst[1]);
    b = b / N;
    double c = pst[0] - a * std::pow(wl[0], 2) - b * wl[0];

    // Wavelength range (400nm -> 1125nm)
    std::array<double, 29> pst_interpolated;
    for (int i = 0; i < 29; i++)
    {
        double wl = 400 + i * 25;
        pst_interpolated[i] = a * std::pow(wl, 2) + b * wl + c;
    }
    return pst_interpolated;
}





std::array<std::vector<std::array<double, 29>>, 5>
StrayLight::getStrayLightAtDetector(
    std::array<std::vector<int>, 5> &rho_a,
    std::array<std::vector<std::array<double, 29>>, 5> &PST,
    std::vector<double> irradiance_alpha)
{

    std::array<std::vector<std::array<double, 29>>, 5> PST_interpolated;

    for (int AZ = 0; AZ < 5; AZ++)
    {
        // AZ = 0
        std::vector<std::array<double, 29>> PST_AZ = PST[0];
        std::vector<int> rho_AZ = rho_a[0];

        std::vector<std::array<double, 29>> PST_interpolated_AZ;
        PST_interpolated_AZ.reserve(PST_AZ.size());

        std::array<std::vector<double>, 29> transposed_PST_interpolated_AZ;

        // Interpolate for every wavelength

        for (int wavelength = 0; wavelength < 29; wavelength++)
        {
            std::vector<std::array<double, 4>> cubicParameters;

            // Do a PCHIP 1-D monotonic cubic interpolation and add the 4
            // parameters (for cubic polynomial) into cubicParameters vector

            double d0 = 0; // Boundary condition of the derivative
                           // at start of interpolation

            for (int angle_idx = 0; angle_idx < rho_AZ.size() - 2; angle_idx++)
            {
                double alpha =
                    double(
                        (rho_AZ[angle_idx + 1] - rho_AZ[angle_idx]) +
                        2 * (rho_AZ[angle_idx + 2] - rho_AZ[angle_idx + 1])) /
                    (3 * (rho_AZ[angle_idx + 2] - rho_AZ[angle_idx]));

                // Boundary condition of the derivative at rho_AZ[angle_idx+1]

                double d1 =
                    (PST_AZ[angle_idx + 1][wavelength] -
                     PST_AZ[angle_idx][wavelength]) *
                    (PST_AZ[angle_idx + 2][wavelength] -
                     PST_AZ[angle_idx + 1][wavelength]) /
                    (alpha *
                         (PST_AZ[angle_idx + 2][wavelength] -
                          PST_AZ[angle_idx + 1][wavelength]) *
                         (rho_AZ[angle_idx + 1] - rho_AZ[angle_idx]) +
                     (1 - alpha) *
                         (PST_AZ[angle_idx + 1][wavelength] -
                          PST_AZ[angle_idx][wavelength]) *
                         (rho_AZ[angle_idx + 2] - rho_AZ[angle_idx + 1]));

                // Get parameters in the region rho_AZ[angle_idx]
                std::array<double, 4> parameters = getCubicParameters(
                    double(rho_AZ[angle_idx]), double(rho_AZ[angle_idx + 1]),
                    PST_AZ[angle_idx][wavelength],
                    PST_AZ[angle_idx + 1][wavelength], d0, d1);
                cubicParameters.push_back(parameters);

                d0 = d1;
            }

            // Extrapolate the PST for every angle in irradiance_alpha
            std::vector<double> extrapolated =
                extrapolate(irradiance_alpha, rho_AZ, cubicParameters);
            transposed_PST_interpolated_AZ[wavelength] = extrapolated;
        }

        // Transpose the transposed_PST_interpolated_AZ

        for (int i = 0; i < transposed_PST_interpolated_AZ[0].size(); i++)
        {
            std::array<double, 29> x;
            for (int j = 0; j < 29; j++)
            {
                x[j] = transposed_PST_interpolated_AZ[j][i];
            }
            PST_interpolated_AZ.push_back(x);
        }

        PST_interpolated[AZ] = PST_interpolated_AZ;
    }

    return PST_interpolated;
}





std::vector<double> StrayLight::extrapolate(std::vector<double> &irradiance_alpha, std::vector<int> &rho, std::vector<std::array<double, 4>> &parameters)
{
    std::vector<double> extrapolated;
    for (double alpha : irradiance_alpha)
    {
        unsigned int idx = 0;
        while (rho[idx+1] < alpha) 
        {
            idx++;
        }

        std::array<double, 4> params = parameters[idx];
        auto f = [params](double a)
        {
            return params[0] * std::pow(a, 3) + params[1] * std::pow(a, 2) +
                   params[2] * a + params[3];
        };

        extrapolated.push_back(f(alpha));
    }

    // return extrapolated
    return extrapolated;
}






std::array<std::vector<double>, 5> StrayLight::getNumberOfStraylightPhotoelectronsAtDetector(std::array<std::vector<std::array<double, 29>>, 5> &PST)
{

    std::array<double, 29> energyOfPhoton;
    std::array<double, 29> wavelengths;
    std::array<std::vector<double>, 5> numberOfStraylightPhotoelectronsAtDetector;

    for (unsigned int wavelength_idx = 0; wavelength_idx < 29; wavelength_idx++)
    {
        double wavelength = 400 + wavelength_idx * 25;
        energyOfPhoton[wavelength_idx] =
            Constants::CLIGHT * Constants::HPLANCK / (wavelength * 1.e-9);
        wavelengths[wavelength_idx] = wavelength;

    }

    for (int i = 0; i < 5; i++)
    {
	std::vector<double> strayLightElectronsPerPixelPerSecond;
        std::vector<std::array<double, 29>> PST_AZ = PST[i];

	for (std::array<double, 29> E_PST : PST_AZ)
        {
            // We integrate out the wavelengths for every gridpoint
	    double integral = 0;
	    for (int wl_idx = 0; wl_idx < 28; wl_idx++)
            {
                integral = integral + 0.5 *
		        (wavelengths[wl_idx + 1] - wavelengths[wl_idx]) *
                        (E_PST[wl_idx] / energyOfPhoton[wl_idx] +
                         E_PST[wl_idx + 1] / energyOfPhoton[wl_idx + 1]);
            }

            double straylightPhotoelectrons =
                integral * std::pow(pixelSize, 2);
            
	    strayLightElectronsPerPixelPerSecond.push_back(straylightPhotoelectrons);
        }
        numberOfStraylightPhotoelectronsAtDetector[i] = strayLightElectronsPerPixelPerSecond;
    }
    return numberOfStraylightPhotoelectronsAtDetector;

}





/**
 * \brief: Returns the 4 parameters that fit the 1-D monotonic cubic
 *         interpolation *in the region x_0 -> x_1.
 *
 * \details: See scipy PchipInterpolator (https://t.ly/8cqwu)
 *
 * \return: Array (parameters) with four parameters that define a cubic
polynomial so that:
 *          f(x) = parameters[0]*x^3 + parameters[1]*x^2 + parameters[2]*x +
parameters[3], and
 *          f(x_0) = y_0; f(x_1) = y_1; f'(x_0) = d0; f'(x_1) = d1
 *
 * \param: pstPath  Path to the corresponding PST file.
 */
std::array<double, 4> StrayLight::getCubicParameters(double x_0, double x_1,
                                                     double y_0, double y_1,
                                                     double d0, double d1)
{

    double a =
        (2 * (y_1 - y_0) - (x_1 - x_0) * (d1 + d0)) / (std::pow(x_0 - x_1, 3));

    double b = ((d1 - d0) + 3 * a * (std::pow(x_0, 2) - std::pow(x_1, 2))) /
               (2 * (x_1 - x_0));
    double c = d0 - 3 * a * std::pow(x_0, 2) - 2 * b * x_0;
    double d = y_1 - a * std::pow(x_1, 3) - b * std::pow(x_1, 2) - c * x_1;

    std::array<double, 4> parameters = {a, b, c, d};

    // std::cout << " " << std::endl;
    // std::cout << a * std::pow(x_0, 3) + b * std::pow(x_0, 2) + c * x_0 + d
    // << " => " << y_0 << std::endl;
    // std::cout << a * std::pow(x_1, 3) + b * std::pow(x_1, 2) + c * x_1 + d << " => " << y_1 << std::endl;
    return parameters;
}





std::tuple<std::vector<double>, std::vector<double>, std::vector<arma::vec>,
            std::vector<arma::vec>>
StrayLight::getIrradianceAtCamera(Camera &camera, std::vector<GridPoint> grid,
                                  std::vector<arma::vec> emmitterIrradiance,
                                  arma::vec emmitterPosition,
                                  arma::vec cameraPosition)
{

    // Get the max angle where light can fall onto the camera

    double cos_alpha_max = cos(radiusFOV);

    // Camera pointing vectors

    double alpha, delta;
    double lambda, beta;

    tie(alpha, delta) = camera.focalPlaneToSkyCoordinates(0, 0);
    equatorial2ecliptic(alpha, delta, lambda, beta);

    

    arma::vec nCamera = {cos(lambda) * cos(beta), sin(lambda) * cos(beta),
                         sin(beta)};
    std::vector<double> irradiance_alpha;
    std::vector<double> grid_alpha;

    std::vector<arma::vec> irradiance_E;
    std::vector<arma::vec> projected_irradiance_E;

    for (unsigned int idx = 0; idx < grid.size(); idx++)
    {

        arma::vec n_Cam_Grid =
            emmitterPosition + grid[idx].point - cameraPosition;
        n_Cam_Grid =
            n_Cam_Grid /
            sqrt(arma::as_scalar(n_Cam_Grid.t() * n_Cam_Grid));

        double cos_alpha1 =
            -1*arma::as_scalar(grid[idx].point.t() * n_Cam_Grid) /
            sqrt(arma::as_scalar(grid[idx].point.t() * grid[idx].point));
        double cos_alpha2 = arma::as_scalar(nCamera.t() * n_Cam_Grid);

        
        arma::vec E(29);

        if (cos_alpha1 < 0)
        {
            E.zeros();

        }
        else if (cos_alpha2 < 0)
        {
            E.zeros();
        }
        else
        {
            E = emmitterIrradiance[idx] * grid[idx].size * cos_alpha1;
        }
	irradiance_E.push_back(E);

        double objectGridCameraAlpha	= 180*acos(cos_alpha1)/Constants::PI;
        double gridCameraPointingAlpha	= 180*acos(cos_alpha2)/Constants::PI;

        arma::vec Eproj(29);
        if (cos_alpha2 < cos_alpha_max)
        {
            Eproj.zeros();
        }
        else
        {
            Eproj = E * cos_alpha2;
        }

        irradiance_alpha.push_back(gridCameraPointingAlpha);
        grid_alpha.push_back(objectGridCameraAlpha);

        projected_irradiance_E.push_back(Eproj);

    }

    return std::tie(irradiance_alpha, grid_alpha, irradiance_E, projected_irradiance_E);
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
