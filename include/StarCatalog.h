#ifndef STARCATALOG_H
#define STARCATALOG_H

#include <string>
#include <vector>
#include <stdexcept>

#include "Logger.h"
#include "Units.h"
#include "SkyCoordinates.h"


using namespace std;



struct StarRecord
{
    public:

        StarRecord(unsigned int starID, double RA, double dec, double Vmag)
        : ID(starID), RA(RA), dec(dec), Vmag(Vmag)
        {};

        ~StarRecord(){};

        StarRecord(StarRecord &&starRecord)
        : ID(starRecord.ID), RA(starRecord.RA), dec(starRecord.dec), Vmag(starRecord.Vmag)
        {};

        const unsigned int   ID;        // Star identification number
        const double         RA;        // Right Ascension [rad]
        const double         dec;       // Declination [rad]
        const double         Vmag;      // Johnson V magnitude
};






class StarCatalog
{
    public:

        StarCatalog();
        StarCatalog(const StarCatalog &starCatalog);
        StarCatalog(StarCatalog &&starCatalog);
        ~StarCatalog();

        long size();
        void addStar(unsigned int starID, double RA, double dec, double Vmag, Unit angleUnit);
        StarRecord operator[](unsigned int index) const;

        StarCatalog getStarsWithinRadiusFrom(double RA0, double dec0, double radius, Unit angleUnit);

    protected:

        long Nstars;
        vector<unsigned int> starID;        // Star identification number
        vector<double> RA;                  // Right Ascension [rad]
        vector<double> dec;                 // Declination [rad]
        vector<double> Vmag;                // Johnson V magnitude

    private:

};



#endif