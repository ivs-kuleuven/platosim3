#include "StrayLight.h"
#include "ConfigurationParameters.h"
#include "Constants.h"
#include "armadillo"
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>
#include <cmath>


#include <chrono>

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
                       Camera &camera, Detector &detector)
  : HDF5Writer(hdf5File), camera(camera), detector(detector)
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

    // Get the coordinates of the telescope reference frame

    arma::vec nx = {1, 0, 0};
    arma::vec ny = {0, 1, 0};
    arma::vec nz = {0, 0, 1};

    telescopeAxis.col(0) = camera.telescopeToSkyCoordinates(nx);
    telescopeAxis.col(1) = camera.telescopeToSkyCoordinates(ny);
    telescopeAxis.col(2) = camera.telescopeToSkyCoordinates(nz);


    double radiusMoon = 1.7381e6; // [m]
    double moon_reflectivity = 0.11;

    // Read in the positions of the file and save them into the vectors.

    std::string orbitPath =
        configParam.getAbsoluteFilename("StrayLight/FilePath");
    readInFile(orbitPath, sc_positions, moon_positions, sun_positions);

    // Read in the PST file and save it into vectors.

    std::string pstPath = configParam.getAbsoluteFilename("StrayLight/PstPath");
    getPST(pstPath);


    // Create the two celestial objects
    moon.radius = radiusMoon;
    moon.reflectivity = moon_reflectivity;
}





/**
 *
 * \brief: Returns the straylight of the moon (in #electrons) in the subfield pixel (row, column)
 *
 * \param: row and column of the subfieldpixel
 */
double StrayLight::getStrayLightMoon(double time)
{

    // We should give the positions of the sc, moon and sun
    arma::vec sun_pos = sun_positions[0];
    arma::vec moon_pos = moon_positions[0];
    arma::vec sc_pos = sc_positions[0];

    
    // Transform these into a reference frame where moon lies in the origin and
    // sun lies on the z-axis.

    sun_pos  = sun_pos - moon_pos;
    sc_pos   = sc_pos- moon_pos;
    moon_pos = moon_pos - moon_pos;

    double A = std::sqrt(sun_pos[0]*sun_pos[0] + sun_pos[1]*sun_pos[1]);
    double N = std::sqrt(sun_pos[0]*sun_pos[0] + sun_pos[1]*sun_pos[1] + sun_pos[2]*sun_pos[2]);

    double x = sun_pos[0];
    double y = sun_pos[1];
    double z = sun_pos[2];
    
    arma::Mat<double> rotation = { { y*N, -x*N, 0},
                                   { x*z,  y*z, -A*A},
                                   { x*A,  y*A,  z*A},};

    rotation = rotation/(A*N);

    sun_pos  = rotation * sun_pos;
    moon_pos = rotation * moon_pos;
    sc_pos   = rotation * sc_pos;


    arma::Mat<double> rotatedAxis = rotation * telescopeAxis;
    return getStrayLightObject(moon, sun_pos, moon_pos, sc_pos, rotatedAxis, 1000);
}






/**
 *
 * \brief: Returns the straylight of the object (in #electrons) in the subfield pixel
 *        (row, column).
 *
 * \param: object        The object from which we want to get the straylight
 * \param: sun_pos       Position of the sun
 * \param: object_pos    Position of the object
 * \param: sc_pos        Position of the spacecraft
 * \param: gridPoints    Amount of gridpoints used to model the object
 * 
 */
double StrayLight::getStrayLightObject(CelestialObject object, arma::vec sun_pos, arma::vec object_pos, arma::vec sc_pos, arma::Mat<double> telescopeAxis, unsigned int nGrid)
{

    std::cout << "Estimaged reflected area: " << 50*(1 - acos(arma::dot(sun_pos,sc_pos) / (arma::norm(sun_pos)*arma::norm(sc_pos))) / 3.1415 ) << "%" << std::endl;
    std::chrono::steady_clock clock;
    auto start_integration = clock.now();

    // This models the total straylight from an object by reflecting light from (nGrid * nGrid) grid points on its surface.
    // Each grid point contributes to the straylight based on three key angles and its surface area.
    // The contribution is calculated as follows:

    // Contribution = Surface Area × cos(gpIrradiance) × cos(gpRadiance) × cos(scIrradiance) [m^2]

    // Where:
    //    -> gpIrradiance: The angle at which sunlight reaches the grid point relative to the normal at that point.
    //    -> gpRadiance: The angle at which the spacecraft is positioned relative to the grid point's normal.
    //    -> scIrradiance: The angle at which reflected light reaches the camera relative to the optical axis.

    // It is important to note that we only consider contributions where the cosine of all angles is positive. 
    // This condition ensures that light can effectively reach the camera, as negative values would indicate
    // angles where light does not contribute to straylight detection.

    double sA = 0; // Reflected exposed surface area
    double tA = 0; // Total exposed surface area
    double maxGridDependencies = 0; // Upper limit to the grid dependencies
    double gridDependencies = 0;    // Total contribution of the grid dependencies
    for (unsigned int i=1; i<=nGrid; i++)
    {
        for (unsigned int j=0; j<=nGrid; j++)
        {
            double theta = i*Constants::PI / (2*nGrid);
            double phi   = 2*j*Constants::PI / nGrid;

            arma::vec n_object = {sin(theta)*cos(phi),
                          sin(theta)*sin(phi),
                          cos(theta)};
            double dA = std::pow(object.radius, 2)*sin(theta)
                        *(Constants::PI / (2*nGrid))* (2*Constants::PI / nGrid);

            double cos_gpIrradiance = arma::dot(sun_pos - object.radius*n_object, n_object) / arma::norm(sun_pos - object.radius*n_object);
            double cos_gpRadiance = arma::dot(sc_pos - object.radius*n_object, n_object) / arma::norm(sc_pos - object.radius*n_object);
            double cos_scIrradiance =
                arma::dot(sc_pos - object.radius * n_object, telescopeAxis.col(2)) /
                arma::norm(sc_pos - object.radius * n_object);

            
            tA += dA;
            if (cos_gpRadiance > 0 && cos_scIrradiance > 0)
            {
                gridDependencies +=
                    cos_gpIrradiance * cos_gpRadiance * cos_scIrradiance *
                    dA;

                sA += dA;
                maxGridDependencies += dA * cos_gpIrradiance;
            }
        }
    }
    
    auto end_integration = clock.now();
    std::chrono::duration<double> elapsed = end_integration-start_integration;
    std::chrono::milliseconds ms_elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(elapsed);
    std::cout << "Exact reflected area: " << (sA / tA)*50 << "%" << std::endl;

    std::cout << "->\tIntegration of gp took: " << ms_elapsed.count() << "ms"
    << std::endl;
    std::cout << "GRID DEPENDENCIES: " << gridDependencies << " m^2\t\tmax: " << maxGridDependencies << " m^2\n" << std::endl;
 

    auto start_pst = clock.now();

    // We now determines the contributions of straylight that are wavelength-dependent. These consist of:
    // - Planck's equation for a black body radiator (B), which models the sunlight that will be reflected. 
    // - The Point Source Transmittance (PST) function, which requires the declination angle (the angle between the optical axis and the incoming ray)
    //   and the azimuth angle (the angle of the projected ray on the camera relative to the x-axis of the telescope). 
    // - Energy per photon, calculated using E = h * c * lambda, where E is energy, h is Planck’s constant, c is the speed of light and lambda is the wavelength.
    
    // The explicit wavelength dependence has already been integrated out (int_lambda B * PST / E). We simply load in the combined effects 
    // from an HDF5 file at the correct declination and azimuth angles. 

    double declination = acos(arma::dot(sc_pos, telescopeAxis.col(2)) / arma::norm(sc_pos));
    declination = rad2deg(declination);

    arma::vec projected = sc_pos - arma::dot(sc_pos, telescopeAxis.col(2))*telescopeAxis.col(2);
    double azimuth = acos(arma::dot(projected, telescopeAxis.col(0)) / arma::norm(projected));
    azimuth = rad2deg(azimuth);

    double pst = getPSTValue(declination, azimuth);
    auto end_pst = clock.now();

    std::chrono::duration<double> elapsed_pst = end_integration-start_integration;
    std::chrono::milliseconds ms_elapsed_pst = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed);
    std::cout << "->\tDetermination of pst took: " << ms_elapsed_pst.count() << "ms" << std::endl;
    std::cout << "PST DEPENDENCIES: " << pst << " ph s-1 m^-2\t\tMax value: " << 1.485e24 << " ph s-1 m^-2"
    << std::endl;
    
    // Lastely we have some constant values that need to be set in order to get
    // the straylight these are:
    // k * pi * pixelSize * exposureTime * Rsun^2 / (d_object_sun^2 *
    // d_object_camera^2), where
    // -> k is the reflexivity of the object
    // -> pixelSize is the area of a pixel (in m)
    // -> exposureTime is the time for one exposure (in s)
    // -> Rsun is the radius of the sun
    // -> d_object_sun is the distance from the object to the sun (in m)
    // -> d_object_camera is the discance from the object to the camera (in m)

    auto object_sun = (sun_pos - object_pos);
    double d_object_sun_sq = arma::dot(object_sun, object_sun);

    auto object_camera = (sc_pos - object_pos);
    double d_object_camera_sq = arma::dot(object_camera, object_camera);


    double const_term = object.reflectivity * std::pow(pixelSize, 2) *
                        std::pow(Constants::SOLARRADIUS, 2) * cycleTime /
                        (d_object_sun_sq * d_object_camera_sq);

    std::cout << "\nCONSTANT DEPENDENCIES: " << const_term << " s" << std::endl;
    return pst*gridDependencies*const_term;


}









double StrayLight::getPSTValue(double declination, double azimuth)
{

    double rho_0 = 0;
    arma::vec pst_azs(4);
    for (int idx_az=0; idx_az < 5; idx_az++)
    {
        for (int idx_rho=0; idx_rho < rhoValues[idx_az].size(); idx_rho++)
        {
            double rho = (rhoValues[idx_az])[idx_rho];

            if (rho < declination)
            {
                rho_0 = rho;
            }
            else
            {
                auto param = (parameters[idx_az])[idx_rho-1];

                auto f = [param](double a)
                {
                    return param[0] * std::pow(a, 3) + param[1] * std::pow(a, 2) +
                           param[2] * a + param[3];
                };

                pst_azs[idx_az] = f(declination);
                break;
            }
      }
    }


    if (azimuth < 0)
    {
        azimuth += 180;
    }
    else if (azimuth > 180)
    {
        azimuth -= 180;
    }

    if (azimuth == 0)
    {
        return pst_azs[0];
    }
    unsigned int idx = 0;
    while (azimuth > azs[idx])
    {
        idx++;
    }

    return ((azimuth - azs[idx-1])*pst_azs[idx] + (azs[idx] - azimuth)*pst_azs[idx-1])/(azs[idx] - azs[idx-1]);
}






/**
 *
 * \brief: Extrapolate the straylight [#e-/s]for the value az. The straylight is
 *         is known for all the values in AZs and given in electronsAtDetector.
 *
 * \param: electronsAtDetector   The values in straylight for the values in AZs.
 * \param: az                    The az value for which we want the straylight.
 *
 * \note: The straylight is known at five points (AZs), we can connect these
 *        points with a fourth order polynomial with coefficients given by
 *        column p. These coefficients are defined such that for any i in
 *        {0,1,2,3,4},
 *        for a[i] = { AZs[i]^4, AZs[i]^3, AZs[i]^2, AZs[i], 1}, we have
 *        a * p = y[i] = electronsAtDetector[i]. This can be expressed as
 *        A*p = y, where the matrix A's ith row is a[i]. Thus p can be found by
 *        p = A^-1 * y. The straylight s at az is then given by:
 *        s = sum_i=0^4 p[i] az^(4-i).
 *       
 */
double StrayLight::getStraylightFromAZ(const std::array<double, 5> &electronsAtDetector, double az)
{

    // We define the matrix A
    arma::mat A(5, 5);
    for (int idx = 0; idx < 5; idx++)
    {
        arma::Row<double> a = { std::pow(AZs[idx], 4), std::pow(AZs[idx], 3), std::pow(AZs[idx], 2), std::pow(AZs[idx], 1), 1};
        A.row(idx) = a;
    }


    // Define the Col y
    arma::Col<double> y = {electronsAtDetector[0], electronsAtDetector[1], electronsAtDetector[2], electronsAtDetector[3], electronsAtDetector[4]};

    // p is given by p = A^-1 * y
    arma::Col<double> p = arma::inv(A) * y;


    arma::Row<double> x = {std::pow(az,4), std::pow(az,3), std::pow(az,2), az, 1};
    return arma::as_scalar(x*p);
}













/**
 * \brief: Calculates the spectral radiance for all gridpoints.
 *
 * \param: sun           Position of the sun.
 * \param: object        Position of the celestial object.
 * \param: reflectivity  Reflectivity of the celestial object.
 * \param: grid          The grid around the celstial object.
 *
 * /note: This method changes the grid variable, so that after the method has
run,
 *        it no longer contains gridpoints that don't receive light from the
sun.
 *
 * /output: objectIrradiance     Irradiance of the celestial object on every
 *                               gridpoint that sees light. [W/m^2*m],
 *                               dim [#lighted_grid_points x #wavelengths]
 *
 */

std::vector<arma::vec> StrayLight::getCelestialObjectGridSpectralRadiance(
    arma::vec sun, arma::vec object, double reflectivity,
    std::vector<GridPoint> &grid)
{

    // Normalized vector from object to sun.
    arma::vec sunObejectVector = (sun - object);
    double distanceSunObject = sqrt(dot(sunObejectVector, sunObejectVector));
    sunObejectVector = sunObejectVector / distanceSunObject;

    // Solar irradiance at the object for different wavelengths
    std::array<double, 29> solarIrradiance =
        solarSpectralIrradiance(distanceSunObject);
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

            // reflectivity of the celestial object
            double a = reflectivity * cos_gamma * cos_gamma / Constants::PI;

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
    std::cout << "Amount of gridpoints: " << grid.size() << " n: " << nPoints << std::endl;
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






void StrayLight::readInFile(std::string orbitPath, std::vector<arma::vec> &sc_position,
                            std::vector<arma::vec> &moon_position,
                            std::vector<arma::vec> &sun_position)
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

            double time = t-t0;
            
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
            sc_position.push_back(sc_row * 1000);
            sun_position.push_back(sun_row * 1000);
            moon_position.push_back(moon_row * 1000);
            times.push_back(time);
        }
    }
}





/**
 * \brief:    Calculates the solar irradiance using Planck's law
 *
 * \param:    distance  distance from object to the sun.
 *
 * \output:   Array containing the irradiance for different wavelengths [W/(m^2*m)]
 *
 */
std::array<double, 29> StrayLight::solarSpectralIrradiance(double distance)
{

    // Wavelength range (400nm -> 1125nm)
    std::array<double, 29> irradiance;

    for (unsigned int i = 0; i < 29; i++)
    {
        double wavelength = (400 + i * 25)*1e-9; // [m]
        // Spectral Solar Radiance
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
 *
 */
void StrayLight::getPST(std::string pstPath)
{
    // Save PST file in datastructure
    if (!FileUtilities::fileExists(pstPath))
    {
        throw FileException("Straylight: trying to load the PST HDF5 file (" + pstPath + "), but file doesn't exist.");
    }
    try
    {
        pstFile.open(pstPath);
    }
    catch (H5::FileIException ex)
    {
        Log.error("H5::FileIException: " + string(ex.getCDetailMsg()));
        throw H5FileException("Straylight: Could not open HDF5 file: " + pstPath);
    }

    std::array<std::vector<double>, 5> pstValues;
    //std::array<std::vector<double>, 5> rhoValues;
    int idx = 0;
    for (int az : azs)
    {
        std::string group = "/az_" + std::to_string(az);
        vector<double> pst;
        vector<double> rho;

        if (pstFile.hasGroup(group))
        {

          // Read in the pst values 

          if (pstFile.hasDataset(group, "pst")) {pstFile.readArray(group, "pst", pst);}
          else
          {
            throw H5FileException("Straylight: PST file does have the dataset pst in group:  " + group);
          }

          // Read in the corresponding rho values

          if (pstFile.hasDataset(group, "rho")) {pstFile.readArray(group, "rho", rho);}
          else
          {
            throw H5FileException("Straylight: PST file does have the dataset rho in group:  " + group);
          }
        }
        pstValues[idx] = pst;
        rhoValues[idx] = rho;
        idx++;
      }
    extrapolatePST(rhoValues, pstValues);
    std::cout << "Done setting up" << std::endl;
}










/**
 * TODO write this
 *
 */
void StrayLight::extrapolatePST(std::array<std::vector<double>, 5> rhoValues,
                                  std::array<std::vector<double>, 5> pstValues)
{

    for (int i = 0; i < 5; i++)
    {
        std::vector<double> PST = pstValues[i];
        std::vector<double> rho = rhoValues[i];

        std::vector<std::array<double, 4>> cubicParameters;

        // Do a PCHIP 1-D monotonic cubic interpolation and add the 4
        // parameters (for cubic polynomial) into cubicParameters vector

        double d0 = 0; // Boundary condition of the derivative
                       // at start of interpolation

        for (int angle_idx = 0; angle_idx < rho.size() - 2; angle_idx++)
        {
            double alpha =
                double(
                       (rho[angle_idx + 1] - rho[angle_idx]) +
                       2 * (rho[angle_idx + 2] - rho[angle_idx + 1])) /
              (3 * (rho[angle_idx + 2] - rho[angle_idx]));

            // Boundary condition of the derivative at rho[angle_idx+1]
            double d1;
            if ((PST[angle_idx] == PST[angle_idx + 1]) &&
                (PST[angle_idx] == PST[angle_idx + 2]))
            {
                d1 = 0;
            }
            else{
                d1 =
                    (PST[angle_idx + 1] -
                     PST[angle_idx]) *
                    (PST[angle_idx + 2] -
                     PST[angle_idx + 1]) /
                    (alpha *
                         (PST[angle_idx + 2] -
                          PST[angle_idx + 1]) *
                         (rho[angle_idx + 1] - rho[angle_idx]) +
                     (1 - alpha) *
                         (PST[angle_idx + 1] -
                          PST[angle_idx]) *
                         (rho[angle_idx + 2] - rho[angle_idx + 1]));
            }
            // Get parameters in the region rho_AZ[angle_idx]
            std::array<double, 4> param = getCubicParameters(
                    double(rho[angle_idx]), double(rho[angle_idx + 1]),
                                            PST[angle_idx], PST[angle_idx + 1], d0, d1);

            cubicParameters.push_back(param);

            d0 = d1;
        }
        parameters[i] = cubicParameters;
    }
}








/**
 *  /brief: extrapolate the values for irradiance_alpha for the piecewise qubic
 *          polynomial defined by parameters defined at the angles rho.
 *
 * \param: irradiance_alpha:  the angles for which we want the obtain the values
 * \param: rho:               the angles that define the intervals where the
 *                            piecwise function is defined.
 * \param: parameters:        the parameters that define the cubic polynoom
 *                            in the intervales given by rho.
 *
 * \note: This function is only used in StrayLight::interpolatePSToverRho
 */
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





/**
 * TODO
 */
std::array<std::vector<double>,5> StrayLight::getNumberOfStraylightPhotoelectronsAtDetector(std::array<std::vector<arma::vec>, 5> &straylight)
{

    arma::vec energyOfPhoton(29); // [J]
    std::array<std::vector<arma::vec>, 5> electronsAtDetectorPerWavelengthPerSecond;   // [#e- / ( s * m)]


    // fill the energyOfPhoton and wavelengths array 
    for (unsigned int wavelength_idx = 0; wavelength_idx < 29; wavelength_idx++)
    {
        double wavelength = 400 + wavelength_idx * 25; // [nm]
        energyOfPhoton[wavelength_idx] =
            Constants::CLIGHT * Constants::HPLANCK / (wavelength * 1.e-9); // [J]
    }

    // get the number of photoelectrons at the detectorfor every wavelength
    for (int az = 0; az < 5; az++)
    {
        for (auto light : straylight[az])
        {
            electronsAtDetectorPerWavelengthPerSecond[az].push_back(
                light / energyOfPhoton);
        }
    }

    std::array<std::vector<double>,5> electronsAtDetectorPerSecond;
    for (int az = 0; az < 5; az++)
    {

	std::vector<double> electronsAtGridpoint;
	for (arma::vec electrons_wl : electronsAtDetectorPerWavelengthPerSecond[az])
        {
            // We integrate out the wavelengths dependency for every gridpoint.
            // We use the composite trapezoidal rule.
	    double integral = 0;
	    for (int wl_idx = 0; wl_idx < 28; wl_idx++)
            {
                integral = integral + 0.5 * (25E-9) *   
		           (electrons_wl[wl_idx] + electrons_wl[wl_idx+1]);
            }
            electronsAtGridpoint.push_back(integral);
        }
        electronsAtDetectorPerSecond[az] = electronsAtGridpoint;
    }
    return electronsAtDetectorPerSecond;
}





/**
 * \brief: Returns the 4 parameters that fit the 1-D monotonic cubic
 *         interpolation *in the region x_0 -> x_1.
 *
 * \details: See scipy PchipInterpolator (https://t.ly/8cqwu)
 *
 * \return: Array (parameters) with four parameters that define a cubic
 *          polynomial so that:
 *          f(x) = parameters[0]*x^3 + parameters[1]*x^2 + parameters[2]*x +
 *          parameters[3], and
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

    auto f = [parameters](double x)
    {
        return  parameters[0] * std::pow(x, 3) + parameters[1] * std::pow(x, 2) +
               parameters[2] * x + parameters[3];
    };
    auto df = [parameters](double x)
    {
        return 3*parameters[0] * std::pow(x, 2) + 2*parameters[1] * x +
               parameters[2];
    };

    return parameters;
}











/**
 * TODO
 */
std::tuple<std::vector<double>, std::vector<double>, std::vector<arma::vec>,
           std::vector<arma::vec>>
StrayLight::getIrradianceAtCamera(Camera &camera, double row, double column,
                                  std::vector<GridPoint> grid,
                                  std::vector<arma::vec> emmitterIrradiance,
                                  arma::vec emmitterPosition,
                                  arma::vec cameraPosition)
{

    // Get the max angle where light can fall onto the camera

    double cos_alpha_max = cos(radiusFOV);

    // pointing vectors
    double xFPmm, yFPmm;
    tie(xFPmm, yFPmm) = detector.pixelToFocalPlaneCoordinates(row, column);
    

    double alpha, delta;
    double lambda, beta;
    tie(alpha, delta) = camera.focalPlaneToSkyCoordinates(xFPmm, yFPmm);
    // tie(alpha, delta) = camera.focalPlaneToSkyCoordinates(0, 0);
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
        double rho_Cam_Grid = sqrt(arma::as_scalar((emmitterPosition + grid[idx].point - cameraPosition).t() * (emmitterPosition + grid[idx].point - cameraPosition)));
        n_Cam_Grid =
            n_Cam_Grid / rho_Cam_Grid;



        double beta0 = asin(((emmitterPosition - cameraPosition)[2]) /
                           sqrt(arma::as_scalar(
                                  (emmitterPosition - cameraPosition).t() *
                                  (emmitterPosition - cameraPosition))));

        double lambda0 = asin(((emmitterPosition - cameraPosition)[1]) /
                         (sqrt(arma::as_scalar(
                              (emmitterPosition - cameraPosition).t() *
                              (emmitterPosition - cameraPosition))) *
                          cos(beta0)));



        double cos_alpha1 =
            -1*arma::as_scalar(grid[idx].point.t() * n_Cam_Grid) /
            sqrt(arma::as_scalar(grid[idx].point.t() * grid[idx].point));
        double cos_alpha2 = arma::as_scalar(nCamera.t() * n_Cam_Grid);

        arma::vec E(29);

        if (cos_alpha1 < 0)
        {
            // Emmitter point doesn't fall on camera
            E.zeros();

        }
        else if (cos_alpha2 < 0)
        {
            // Emmitted light doesn't fall on camera
            E.zeros();
        }
        else
        {
            E = emmitterIrradiance[idx] * grid[idx].size * cos_alpha1 /
                (rho_Cam_Grid * rho_Cam_Grid);
        }

	irradiance_E.push_back(E);

        double gridCameraAlpha	= 180*acos(cos_alpha1)/Constants::PI;
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
        grid_alpha.push_back(gridCameraAlpha);
	
        projected_irradiance_E.push_back(Eproj);

    }

    return std::tie(irradiance_alpha, grid_alpha, irradiance_E, projected_irradiance_E);
}




/**
 * TODO
 */
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
