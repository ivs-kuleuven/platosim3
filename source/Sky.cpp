
#include "Sky.h"




/**
 * \brief Constructor
 * 
 * \param configParams    Configuration parameters as read from the (e.g. yaml) inputfile
 * 
 */

Sky::Sky(ConfigurationParameters &configParams)
{
    // Configure this Sky 

    configure(configParams);

    // Open and read the file containing the position and magnitude of all stars
    // The path of the starInputfile should have been set in configure().

    ifstream myfile(starInputfile);
    if (myfile.is_open())
    {
        string temp;
        unsigned int n = 0;
        while (getline(myfile, temp))
        {
            istringstream buffer(temp);
            vector<double> numbers((istream_iterator<double>(buffer)), istream_iterator<double>());
            starCatalog.addStar(n, numbers[0], numbers[1], numbers[2], Angle::degrees);
            n++;
        }

        myfile.close();

        Log.info("Sky: found " + to_string(starCatalog.size()) + " stars in file " + starInputfile);
    }
    else
    {
        Log.error("Sky: Cannot open star catalog file " + starInputfile);
        exit(1);
    }
}








/**
 * \brief Destructor
 */

Sky::~Sky()
{

}









/**
 * \brief Configure the Sky class with the user given input parameters
 */

void Sky::configure(ConfigurationParameters &configParams)
{
    starInputfile = configParams.getAbsoluteFilename("ObservingParameters/StarCatalogFile");
}










/**
 * \brief  Return the equatorial sky coordinates of the star with a given ID.
 * 
 * \param id               Sequential number of the star 
 * \param outputAngleUnit  Either Angle::degrees or Angle::radians
 * 
 * \return (RA, dec)  Equatorial coordinates of the star [rad]
 */

pair<double, double> Sky::getCoordinatesOfStarWithID(int id, Unit outputAngleUnit)
{
    if ((id < 0) || (id >= starCatalog.size()))
    {
        string errorMessage = "Sky: getStarCoordinatesOfStarWithID(): id " + to_string(id) 
                            + " is not in {0,.., " + to_string(starCatalog.size()-1) + "}";
        Log.error(errorMessage);
        throw IllegalArgumentException("errorMessage");
    }
    else
    {
        const auto star = starCatalog[id];
        return make_pair(star.RA * outputAngleUnit, star.dec * outputAngleUnit);
    }
}










/**
 * \brief  Return the V magnitude of the star with given ID
 * 
 * \param id       Sequential number of the star 
 * \return Vmag    V-magnitude of the star
 */

double Sky::getVmagnitudeOfStarWithID(int id)
{
    if ((id < 0) || (id >= starCatalog.size()))
    {
        string errorMessage = "Sky: getVmagnitudeOfStarWithID(): id " + to_string(id) 
                            + " is not in {0,.., " + to_string(starCatalog.size()-1) + "}";
        Log.error(errorMessage);
        throw IllegalArgumentException("errorMessage");
    }
    else
    {
        const auto star = starCatalog[id];
        return star.Vmag;
    }
}












/**
 * \brief  Given a circle on the sky, return a catalog with all stars from the database within that circle.
 * 
 * \note The stars right on the circle are also included in the catalog.
 * 
 * \param RA0        Right Ascencsion of the center point of the circle on the sky
 * \param dec0       Declination of the center point of the circle on the sky
 * \param radius     Radius of the circle on the sky
 * \param angleUnit  If the input angles are in degrees: Angle:degrees, if in radians: Angle::radians 
 * 
 * \return           A StarCatalog object containing the ID, RA, dec, and Vmag of each star within the circle.
 */

StarCatalog Sky::getStarsWithinRadiusFrom(double RA0, double dec0, double radius, Unit angleUnit)
{
    return starCatalog.getStarsWithinRadiusFrom(RA0, dec0, radius, angleUnit);
}











/**
 * \brief Return the solar radiant flux at the given wavelength, measured above the atmosphere of the earth.
 * 
 * \details This function uses the Wehrli (1985) solar irradiance table
 *          plus linear interpolation. Note that the units of the tabulated
 *          data radiant fluxes are \f$ J s^{-1} m^{-2} (nm)^{-1} \f$, where nm
 *          is nanometer (unit of wavelength) while the units of the radiant
 *          flux that this function returns is SI: \f$ J s^{-1} m^{-2} m^{-1} \f$
 * 
 * \param lambda  Wavelength [m], should be in [199.5 nm, 100075 nm]
 * 
 * \return solar radiant flux at air mass zero [\f$J s^{-1} m^{-2} m^{-1}\f$]
 */

double Sky::solarRadiantFlux(double lambda)
{
    int i;         // Index of the first point, defining the linear relation
    double result; // solar radiant flux in SI units.

    // The tabulated wavelengths are in nm, so we have to convert from [m] to [nm]

    double lam = lambda * 1.0e+9;

    // Find the location of the wavelength point lam in the wavelength
    // array skydata::sunwavel[]

    locate(lam, skydata::sunwavel, 920, i);

    // Check if the location was found

    if (i == -1)
    {
        Log.error("Sky::solarRadiantFlux()): no data for the given wavelength");
        exit(1);
    }

    // Do the linear interpolation

    result = skydata::sunflux[i]
                + (skydata::sunflux[i+1] - skydata::sunflux[i])
                / (skydata::sunwavel[i+1] - skydata::sunwavel[i])
                * (lam - skydata::sunwavel[i]);

    // The tabulated flux data are per nm (unit of wavelength), so we should
    // convert to flux per m.

    return result * 1.0e9;             // conversion: (nm)^{-1} -> m^{-1}
}














/**
 * \brief Computes the solar radiant flux between the wavelengths lambda1 and lambda2.
 * 
 * \details This function uses the Wehrli (1985) solar irradiance table
 *          plus the SolarRadiantFlux(lambda) function.
 *          . lambda1 and lambda2 should be between 199.5 nm and 10075 nm.
 *          . This function is overloaded.
 * 
 * \param lambda1  Lower wavelength [m] of the interval, should be in [199.5 nm, 100075 nm]
 * \param lambda2  Upper wavelength [m] of the interval, should be in [199.5 nm, 100075 nm]
 * 
 * \return Integrated solar radiant flux [\f$J s^{-1} m^{-2}\f$]
 */

double Sky::solarRadiantFlux(double lambda1, double lambda2)
{
    double lam1;
    double lam2;

    // Check if the wavelengths are within the table boundaries.

    if ((lambda1 >= 199.5e-9) && (lambda1 <= 10075.0e-9) && (lambda2 >= 199.5e-9) && (lambda2 <= 10075.0e-9))
    {
        // If the integral boundaries are equal, the integral is zero

        if (lambda1 == lambda2)
        {
            return (0.0);
        }

        // Check if the first wavelength is indeed greater than the second

        if (lambda1 < lambda2)
        {
            lam1 = lambda1;
            lam2 = lambda2;
        }
        else
        {
            lam1 = lambda2;
            lam2 = lambda1;
        }

        // Integrate with the trapezium method (Numerical Recipes, pg. 137)

        const int JMAX = 30;
        const double EPS = 1.0e-5;
        double x, tnm, sum, del, olds;
        double s = 0.0;
        int it, j, k;

        olds = -1.0e30;

        for (j = 1; j <= JMAX; j++)
        {
            if (j == 1)
            {
                s = 0.5 * (lam2 - lam1) * (solarRadiantFlux(lam1) + solarRadiantFlux(lam2));
            }
            else
            {
                for (it = 1, k = 1; k < j - 1; k++) it <<= 1;
                tnm = it;
                del = (lam2 - lam1) / tnm;
                x = lam1 + 0.5 * del;
                for (sum = 0.0, k = 1; k <= it; k++, x += del)
                {
                    sum += solarRadiantFlux(x);
                }
                s = 0.5 * (s + (lam2 - lam1) * sum / tnm);
            }

            if (j > 5)
            {
                if (fabs(s - olds) < EPS * fabs(olds) || (s == 0.0 && olds == 0.0)) 
                return s;
            }

            olds = s;
        }

        Log.error("Sky::solarRadiantFlux(): Integration not converged");
        exit (1);
    }
    else
    {
        // This is the case that the given wavelengths were not between the table boundaries.

        Log.error("Sky::solarRadiantFlux(): wavelength must be in [199.5e-9, 10075.0e-9]");
        exit (1);
    }
} 











/**
 * \brief Compute the Solar radiant flux in the given passband
 * 
 * \details This function uses the Wehrli (1985) solar irradiance table
 *          plus the SolarRadiantFlux(lambda) function.
 * 
 * \param lambda     Wavelengths of the passband [m], should be in [199.5 nm, 100075 nm]
 * \param throughput Throughput of the passband
 * 
 * \return Solar radiant flux  [\f$J s^{-1} m^{-2}\f$]
 */

double Sky::solarRadiantFlux(vector<double> &lambda, vector<double> &throughput)
{
    const double lambda1 = lambda[0];
    const double lambda2 = lambda[lambda.size()-1];

    // Check if the passband wavelengths are within [199.5, 10075] nm.

    if ((lambda1 < 199.5e-9) || (lambda1 > 10075.0e-9) || (lambda2 < 199.5e-9) || (lambda2 > 10075.0e-9))
    {
        Log.error("Sky::solarRadiantFlux(): Passband wavelengths not in [199.5, 10075] nm.");
        exit(1);
    }

  // Build up the function you want to integrate

   integrand.clear();
   integrand.resize(lambda.size());

   for (unsigned int i = 0; i < lambda.size(); i++)
   {
      integrand[i] = solarRadiantFlux(lambda[i]) * throughput[i];
   }

   tabfunction.init(lambda, integrand, lambda.size());
   tabfunction.setInterpolationMethod(Linear_Interpolation);

   // Integrate over the throughput band

   const double integral = tabfunction.integrate(lambda[0], lambda[lambda.size()-1]);

   return integral;

}










/**
 * \brief Compute the zodiacal background flux in the wavelength interval [lambda1, lambda2], for a given position in the sky.
 * 
 * \note Some parts of the sky cannot be sampled!
 * 
 * \param alpha    Right ascension coordinate  [rad]
 * \param delta    Declination coordinate      [rad]
 * \param lambda1  Lower wavelength of the interval [m]
 * \param lambda2  Upper wavelength of the interval [m]
 * 
 * \return Zodiacal flux [\f$J s^{-1} m^{-2} sr^{-1}\f$]
 */

double Sky::zodiacalFlux(double alpha, double delta, double lambda1, double lambda2)
{
    double lam, beta;
    double flux500;
    int lam_index, beta_index;

    // Convert the equatorial coordinates to geocentric ecliptic coordinates.
    // All coordinates are in radians.

    auto skyPoint = SkyCoordinates(alpha, delta, Angle::radians);
    tie(lam, beta) = skyPoint.toEcliptic(Angle::radians);

    // The zodiacal light is approximately symmetric with respect to the
    // ecliptic, and with respect to the helioecliptic meridian
    // (= sun-ecliptic poles-antisolar point).

    beta = fabs (beta);
    if (lam > Constants::PI) lam = 2.0 * Constants::PI - lam;

    // Convert from radians to degrees

    beta = rad2deg(beta);
    lam  = rad2deg(lam);

    // Locate the coordinates lam & beta in the coordinate arrays
    // Check if the coordinates are out of boundary.

    locate(lam, skydata::zodlong, 19, lam_index);
    locate(beta,skydata::zodlat,  10, beta_index);

    if ((lam_index == -1) || (beta_index == -1))
    {
        string position = "(" + to_string(rad2deg(alpha)) + ", " + to_string(rad2deg(delta)) + ")";
        Log.warning("Sky::zodiacalFlux(): No data for (alpha, delta) = " + position);
        return 0.0;
    }

    // Check if we don't happen to be in a "hole" in the table

    if (  (skydata::zod[lam_index][beta_index] == -1) || (skydata::zod[lam_index][beta_index+1] == -1)
        || (skydata::zod[lam_index+1][beta_index] == -1) || (skydata::zod[lam_index+1][beta_index+1] == -1))
    {
        string position = "(" + to_string(rad2deg(alpha)) + ", " + to_string(rad2deg(delta)) + ")";
        Log.warning("Sky::zodiacalFlux(): No data for (alpha, delta) = " + position);
        return 0.0;
    }

    // Do a bilinear interpolation

    double dx = (skydata::zodlat[beta_index+1] - beta) / (skydata::zodlat[beta_index+1] - skydata::zodlat[beta_index]);
    double dy = (skydata::zodlong[lam_index+1] - lam) / (skydata::zodlong[lam_index+1] - skydata::zodlong[lam_index]);
    double P1 = skydata::zod[lam_index][beta_index];
    double P2 = skydata::zod[lam_index][beta_index+1];
    double P3 = skydata::zod[lam_index+1][beta_index];
    double P4 = skydata::zod[lam_index+1][beta_index+1];

    flux500 =  dx * dy * P1 + (1.0 - dx) * dy * P2 + dx * (1.0 - dy) * P3  + (1.0 - dx) * (1.0 - dy) * P4;

    // The tabulated fluxes are given in:
    //    10^{-8} J s^{-1} m^{-2} sr^{-1} (micrometer)^{-1}
    // so we need to convert to:
    //    J s^{-1} m^{-2} sr^{-1} m^{-1}

    flux500 *= 0.01;

    // Now we have the zodiacal light at 500 nm, from which we need to
    // derive the zodiacal light flux in the interval [lambda1, lambda2].
    // For this we use the fact that the zodiacal flux has a solar
    // wavelength dependence.

    return (flux500 * solarRadiantFlux(lambda1, lambda2) / solarRadiantFlux(500e-9));
}













/**
 * \brief Compute the zodiacal background flux in the given passband, for a given position in the sky.
 * 
 * \note Some parts of the sky cannot be sampled!
 * 
 * \param alpha       Right ascension coordinate [rad]
 * \param delta       Declination coordinate [rad]
 * \param lambda      Wavelengths of the passband [m]
 * \param throughput  Throughput of the passband
 * 
 * \return  Zodiacal flux [\f$J s^{-1} m^{-2} sr^{-1}\f$]
 */

double Sky::zodiacalFlux(double alpha, double delta, vector<double> &lambda, vector<double> &throughput)
{
    double lam, beta;
    double flux500;
    int lam_index, beta_index;


    // Convert the equatorial coordinates to geocentric ecliptic coordinates.
    // All coordinates are in radians.

    auto skyPoint = SkyCoordinates(alpha, delta, Angle::radians);
    tie(lam, beta) = skyPoint.toEcliptic(Angle::radians);

    // The zodiacal light is approximately symmetric with respect to the
    // ecliptic, and with respect to the helioecliptic meridian (= sun-ecliptic poles-antisolar point).

    beta = fabs(beta);
    if (lam > Constants::PI) lam = 2.0 * Constants::PI - lam;

    // Convert from radians to degrees

    beta = rad2deg(beta);
    lam  = rad2deg(lam);

    // Locate the coordinates lam & beta in the coordinate arrays
    // Check if the coordinates are out of boundary.

    locate(lam,  skydata::zodlong, 19, lam_index);
    locate(beta, skydata::zodlat,  10, beta_index);

    if ((lam_index == -1) || (beta_index == -1))
    {
        string position = "(" + to_string(rad2deg(alpha)) + ", " + to_string(rad2deg(delta)) + ")";
        Log.warning("Sky::zodiacalFlux(): No data for (alpha, delta) = " + position);
        return 0.0;
    }

    // Check if we don't happen to be in a "hole" in the table

    if (  (skydata::zod[lam_index][beta_index] == -1) || (skydata::zod[lam_index][beta_index+1] == -1)
        || (skydata::zod[lam_index+1][beta_index] == -1) || (skydata::zod[lam_index+1][beta_index+1] == -1))
    {
        string position = "(" + to_string(rad2deg(alpha)) + ", " + to_string(rad2deg(delta)) + ")";
        Log.warning("Sky::zodiacalFlux(): No data for (alpha, delta) = " + position);
        return 0.0;
    }

    // Do a bilinear interpolation

    double dx = (skydata::zodlat[beta_index+1] - beta) / (skydata::zodlat[beta_index+1] - skydata::zodlat[beta_index]);
    double dy = (skydata::zodlong[lam_index+1] - lam) / (skydata::zodlong[lam_index+1] - skydata::zodlong[lam_index]);
    double P1 = skydata::zod[lam_index][beta_index];
    double P2 = skydata::zod[lam_index][beta_index+1];
    double P3 = skydata::zod[lam_index+1][beta_index];
    double P4 = skydata::zod[lam_index+1][beta_index+1];

    flux500 =  dx * dy * P1 + (1.0 - dx) * dy * P2 + dx * (1.0 - dy) * P3  + (1.0 - dx) * (1.0 - dy) * P4;

    // The tabulated fluxes are given in:
    //    10^{-8} J s^{-1} m^{-2} sr^{-1} (micrometer)^{-1}
    // so we need to convert to:
    //    J s^{-1} m^{-2} sr^{-1} m^{-1}

    flux500 *= 0.01;

    // Now we have the zodiacal light at 500 nm, from which we need to
    // derive the zodiacal light flux in the passband 'throughput'.
    // For this we use the fact that the zodiacal flux has a solar
    // wavelength dependence.

    return (flux500 * solarRadiantFlux(lambda, throughput) / solarRadiantFlux(500e-9));

}












/**
 * \brief Return the Stellar background (unresolved stars + diffuse galactic background + extragalactic background) 
 *        for the given equatorial coordinates, in the wavelength interval [lambda1, lambda2].
 * 
 * \details Tabulated values of the Pioneer 10 blue passband (spectral range: [395, 495] nm), 
 *          and the Pioneer 10 red passband (spectral range: [590, 690] nm) are used. The value 
 *          in the given interval will be computer by (crude!) linear inter- and extrapolation which.
 *          Care is taken that in the Pioneer 10 passbands, the interpolation returns exactly the 
 *          tabulated values. 
 * 
 * \note Some parts of the sky cannot be sampled!
 * 
 * \param RA       Equatorial coordinate: right ascension [rad]
 * \param dec      Equatorial coordinate: declination [rad]
 * \param lambda1  Begin wavelength of the interval [m]
 * \param lambda2  End   wavelength of the interval [m]
 * 
 * \return Stellar background flux in the Pioneer 10 blue/red passband [\f$J s^{-1} m^{-2} sr^{-1}\f$]
 */

double Sky::stellarBackgroundFlux(double RA, double dec, double lambda1, double lambda2)
{
    double alpha, delta;
    int alpha_index, delta_index;
    double blueflux, redflux;
    double a, b;

    // Convert from radians to degrees

    alpha = rad2deg(RA);
    delta = rad2deg(dec);

    // Locate the coordinates alpha & delta in the coordinate arrays
    // Check if the coordinates are out of boundary.

    locate(alpha, skydata::skyRA, 37, alpha_index);
    locate(delta, skydata::skydec, 25, delta_index);

    if ((alpha_index == -1) || (delta_index == -1))
    {
        string skyPosition = "(" + to_string(alpha) + ", " + to_string(delta) + ") deg";
        Log.warning("Sky::stellarBgFlux(): No data for (alpha, delta) = " + skyPosition);
        return 0.0;
    }


    // Check if we don't happen to be in a "hole" in the table
    // These "holes" are the same for the blue and red passband.

    if (   (skydata::skyblue[alpha_index][delta_index] == -1) || (skydata::skyblue[alpha_index][delta_index+1] == -1)
        || (skydata::skyblue[alpha_index+1][delta_index] == -1) || (skydata::skyblue[alpha_index+1][delta_index] == -1))
    {
        string skyPosition = "(" + to_string(alpha) + ", " + to_string(delta) + ") deg";
        Log.warning("Sky::stellarBgFlux(): No data for (alpha, delta) = " + skyPosition);
        return 0.0;
    }


    // Do a bilinear interpolation for both the blue and the red passband

    double dx = (skydata::skydec[delta_index+1] - delta) / (skydata::skydec[delta_index+1] - skydata::skydec[delta_index]);
    double dy = (skydata::skyRA[alpha_index+1] - alpha) / (skydata::skyRA[alpha_index+1] - skydata::skyRA[alpha_index]);

    double P1 = skydata::skyblue[alpha_index][delta_index];
    double P2 = skydata::skyblue[alpha_index][delta_index+1];
    double P3 = skydata::skyblue[alpha_index+1][delta_index];
    double P4 = skydata::skyblue[alpha_index+1][delta_index+1];

    blueflux =   dx * dy * P1 + (1.0 - dx) * dy * P2 + dx * (1.0 - dy) * P3  + (1.0 - dx) * (1.0 - dy) * P4;

    P1 = skydata::skyred[alpha_index][delta_index];
    P2 = skydata::skyred[alpha_index][delta_index+1];
    P3 = skydata::skyred[alpha_index+1][delta_index];
    P4 = skydata::skyred[alpha_index+1][delta_index+1];

    redflux =   dx * dy * P1 + (1.0 - dx) * dy * P2 + dx * (1.0 - dy) * P3  + (1.0 - dx) * (1.0 - dy) * P4;


    // Convert the result from S10sun units (brightness of 10th magnitude
    // solar type stars per degree square) to SI units:
    //   J s^{-1} m^{-2} sr^{-1}

    blueflux *= 1.2084e-9;
    redflux  *= 1.0757e-9;


    // Do an extra/interpolation to the user given wavelength band
    // Note: Our 'a' is in the formulae 'a/2'.

    const double R1 = 590.0e-9;   // Lower wavelength of Pioneer 10 Red Passb.
    const double R2 = 690.0e-9;   // Upper wavelength of Pioneer 10 Red Passb.
    const double B1 = 395.0e-9;   // Lower wavelength of Pioneer 10 Blue Passb.
    const double B2 = 495.0e-9;   // Upper wavelength of Pioneer 10 Blue Passb.

    a =  (redflux * (B2 - B1) - blueflux * (R2 - R1))
            / ((R2 * R2 - R1 * R1) * (B2 - B1) - (B2 * B2 - B1 * B1) * (R2 - R1));

    b =  (blueflux * (R2 * R2 - R1 * R1) - redflux * (B2 * B2 - B1 * B1))
            / ((R2 * R2 - R1 * R1) * (B2 - B1) - (B2 * B2 - B1 * B1) * (R2 - R1));


    return (a * (lambda2*lambda2 - lambda1*lambda1) + b * (lambda2 - lambda1));
}














/**
 * \brief Approximate the stellar background (unresolved stars + diffuse galactic background 
 *       + extragalactic background) for the given equatorial coordinates, in the given passband.
 *       
 * \details Tabulated values for the Pioneer 10 blue passband (spectral range: [395, 495] nm), 
 *          and for the Pioneer 10 red passband (spectral range: [590, 690] nm) are used.
 *          For lambda <= 690 nm, first the monochromatic flux function is approximated as a 
 *          linear function (increasing with wavelength because there is more red than blue 
 *          background light). Care is taken that the approximated flux is never negative. 
 *          For lambda >= 690 nm, the monochromatic flux is approximated as a constant function 
 *          with the value of the flux at lambda = 690 nm. Care is taken that in the Pioneer 10 
 *          passbands, the interpolation returns exactly the tabulated values.
 *          
 * \note Some parts of the sky cannot be sampled!
 *   
 * \param RA          Equatorial coordinate: right ascension [radians]
 * \param dec         Equatorial coordinate: declination [radians]
 * \param lambda      Wavelength values of the passband  [m]
 * \param throughput  Throughput of the passband
 * 
 * \return  Stellar background flux [\f$J s^{-1} m^{-2} sr^{-1}\f$]
 */

double Sky::stellarBackgroundFlux (double RA, double dec, vector<double> &lambda, vector<double> &throughput)
{
    int alpha_index, delta_index;
    double blueflux, redflux;
    double a, b;

    // Convert from radians to degrees

    double alpha = rad2deg(RA);
    double delta = rad2deg(dec);

    // Locate the coordinates alpha & delta in the coordinate arrays
    // Check if the coordinates are out of boundary.

    locate(alpha, skydata::skyRA, 37, alpha_index);
    locate(delta, skydata::skydec, 25, delta_index);

    if ((alpha_index == -1) || (delta_index == -1))
    {
        string skyPosition = "(" + to_string(alpha) + ", " + to_string(delta) + ") deg";
        Log.warning("Sky::stellarBgFlux(): No data for (alpha, delta) = " + skyPosition);
        return 0.0;
    }


    // Check if we don't happen to be in a "hole" in the table
    // These "holes" are the same for the blue and red passband.

    if (  (skydata::skyblue[alpha_index][delta_index] == -1) || (skydata::skyblue[alpha_index][delta_index+1] == -1)
        || (skydata::skyblue[alpha_index+1][delta_index] == -1) || (skydata::skyblue[alpha_index+1][delta_index] == -1))
    {
        string skyPosition = "(" + to_string(alpha) + ", " + to_string(delta) + ") deg";
        Log.warning("Sky::stellarBgFlux(): No data for (alpha, delta) = " + skyPosition);
        return 0.0;
    }


    // Do a bilinear interpolation for both the blue and the red passband

    double dx = (skydata::skydec[delta_index+1] - delta) / (skydata::skydec[delta_index+1] - skydata::skydec[delta_index]);
    double dy = (skydata::skyRA[alpha_index+1] - alpha) / (skydata::skyRA[alpha_index+1] - skydata::skyRA[alpha_index]);

    double P1 = skydata::skyblue[alpha_index][delta_index];
    double P2 = skydata::skyblue[alpha_index][delta_index+1];
    double P3 = skydata::skyblue[alpha_index+1][delta_index];
    double P4 = skydata::skyblue[alpha_index+1][delta_index+1];

    blueflux =   dx * dy * P1 + (1.0 - dx) * dy * P2 + dx * (1.0 - dy) * P3  + (1.0 - dx) * (1.0 - dy) * P4;

    P1 = skydata::skyred[alpha_index][delta_index];
    P2 = skydata::skyred[alpha_index][delta_index+1];
    P3 = skydata::skyred[alpha_index+1][delta_index];
    P4 = skydata::skyred[alpha_index+1][delta_index+1];

    redflux =   dx * dy * P1 + (1.0 - dx) * dy * P2 + dx * (1.0 - dy) * P3  + (1.0 - dx) * (1.0 - dy) * P4;


    // Convert the result from S10sun units (brightness of 10th magnitude
    // solar type stars per degree square) to SI units:
    //   J s^{-1} m^{-2} sr^{-1}

    blueflux *= 1.2084e-9;
    redflux  *= 1.0757e-9;

    // Compute the linear relation f(lambda) = a * lambda + b
    // to be used for lambda <= 690 nm.

    const double R1 = 590.0e-9;   // Lower wavelength of Pioneer 10 Red Passb.
    const double R2 = 690.0e-9;   // Upper wavelength of Pioneer 10 Red Passb.
    const double B1 = 395.0e-9;   // Lower wavelength of Pioneer 10 Blue Passb.
    const double B2 = 495.0e-9;   // Upper wavelength of Pioneer 10 Blue Passb.

    a =  2.0 * (redflux * (B2 - B1) - blueflux * (R2 - R1))
        / ((R2 * R2 - R1 * R1) * (B2 - B1) - (B2 * B2 - B1 * B1) * (R2 - R1));

    b =  (blueflux * (R2 * R2 - R1 * R1) - redflux * (B2 * B2 - B1 * B1))
        / ((R2 * R2 - R1 * R1) * (B2 - B1) - (B2 * B2 - B1 * B1) * (R2 - R1));


    // Check if the monochromatic flux is indeed increasing for lambda <= 690 nm.

    if (a <= 0.0)
    {
        Log.warning("Sky::stellarBgFlux(): Bad behaviour of monochromatic background flux for lambda <= 690 nm.");
    }

    // Compute where this linear function is zero. Useful to know, because
    // we don't want negative sky background fluxes.

    double root = - b / a;

    // Set up the integrand

    integrand.clear();
    integrand.resize (lambda.size());

    for (unsigned int i = 0; i < lambda.size(); i++)
    {
        if (lambda[i] <= root)
        {
            integrand[i] = 0.0;
        }

        if ((lambda[i] > root) && (lambda[i] <= R2))
        {
            integrand[i] = (a * lambda[i] + b) * throughput[i];
        }

        if (lambda[i] > R2)
        {
            integrand[i] = (a * R2 + b) * throughput[i];
        }
    }

    // Return the integral

    tabfunction.init(lambda, integrand, lambda.size());
    tabfunction.setInterpolationMethod(Linear_Interpolation);
   
    return (tabfunction.integrate(lambda[0], lambda[lambda.size()-1]));
}













/**
 * \brief Given an array[0..N-1] of ascending values, and a value x,
 *        return an index so that array[index] <= x <= array[index+1]
 *        If the value x is out of the array boundaries, index will be set to -1.
 * 
 * \note NO error trapping!
 * 
 * \param x
 * \param array
 * \param N
 * \param index
 */

void Sky::locate(double x, const double *array, int N, int &index)
{
    int index1, index2;

    // Check if x is out of the tabulated range

    if ((x < array[0]) || (x > array[N-1]))
    {
        index = -1;
        return;
    }

    // Check if x happens to be the first element of the array

    if (x == array[0])
    {
        index = 0;
    }

    // Check if x happens to be the last element of the array

    if (x == array[N-1])
    {
        index = N - 2;
    }

    // Find the location with bisection

    index1 = 0;                       // We already checked the lower and upper
    index2 = N - 1;                   // borders.

    unsigned int middle;               // Middle point

    while (index2 - index1 > 1)
    {
        middle = (index2 + index1) >> 1;

        if (x >= array[middle])
        {
            index1 = middle;
        }
        else
        {
            index2 = middle;
        }
    }

    index = index1;
}


