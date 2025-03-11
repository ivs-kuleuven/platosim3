
- File naming: in English its' called "Straylight". The name StrayLight with capital L is therefore confusing. Better: "Straylight.cpp" and
    "Straylight.h"

- The code has a rather strange architecture, which makes it more challenging to debug than it needs to be.

- Much of the code can move from StrayLight.* to a StraylightObject.



class StraylightObject {

    public:

        StraylightObject(double radius, double reflectivity, unsigned int gridSize, string::orbitFilePath);
        irradiance reflectionTowardsSpacecraft(time, position_spacecraft, position_sun);

    protected:
        
        void readOrbitFile(std::string orbitFilePath);

    private:

        double radius;
        double reflectivity;
        arma::vec gridPoints;
        double gridPointArea;
        vec<double> time;
        vec<double> orbit_x;            // [m]
        vec<double> orbit_y;            // [m]
        vec<double> orbit_z;            // [m]
};





- Can we approximate the integral over the moon's surface and/or the integral
  over the wavelengths?

- A StraylightObject could be configured in a Camera object. 

- Grid interpolation and extrapolation should be done in a separate file,
    because this may also be useful in other parts of the PlatoSim code.



In the input yaml file:

Straylight:
    CelestialObjects:
        Moon:
            radius: xxx
            reflectivity: xxx
            orbitFilePath: zzz
        Earth:
            radius: xxx
            reflectivity: yyy
            orbitFilePath: zzz

Camera;
    PstFilePath: inputfiles/Plato_Straylight_PST.txt   # Point Source Transmittance


