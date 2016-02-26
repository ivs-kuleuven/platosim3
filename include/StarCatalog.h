#ifndef STARCATALOG_H
#define STARCATALOG_H

#include <string>
#include <vector>

#include "Logger.h"
#include "Units.h"
#include "SkyCoordinates.h"


using namespace std;



struct StarRecord
{
    public:

        StarRecord(const long starID, const double RA, const double dec, const double Vmag)
        : ID(starID), RA(RA), dec(dec), Vmag(Vmag)
        {};

        ~StarRecord(){};

        StarRecord(StarRecord &&starRecord)
        : ID(starRecord.ID), RA(starRecord.RA), dec(starRecord.dec), Vmag(starRecord.Vmag)
        {};

        const long   ID;        // Star identification number
        const double RA;        // Right Ascension [rad]
        const double dec;       // Declination [rad]
        const double Vmag;      // Johnson V magnitude
};






class StarCatalog
{
    public:

        StarCatalog();
        StarCatalog(const StarCatalog &starCatalog);
        StarCatalog(StarCatalog &&starCatalog);
        ~StarCatalog();

        long size();
        void addStar(const long starID, const double RA, const double dec, const double Vmag, Unit angleUnit);
        StarRecord operator[](long index) const;

        StarCatalog getStarsWithinRadiusFrom(const double RA0, const double dec0, const double radius, Unit angleUnit);

    protected:

        long Nstars;
        vector<long> starID;        // Star identification number
        vector<double> RA;          // Right Ascension [rad]
        vector<double> dec;         // Declination [rad]
        vector<double> Vmag;        // Johnson V magnitude

    private:

};



#endif