
#ifndef DETECTOR_H
#define DETECTOR_H

#include <string>
#include "telescope.h"
#include "subfield.h"


using namespace std;


class Detector
{
    public:

        Detector(ConfigurationParameters configurationParameters);
        ~Detector();

        virtual void takeExposure();
        SubField getSubField();

    protected:

        virtual void integrateLight();
        virtual void readOut();
        virtual void applyOpenShutterSmearing();
        virtual void applyFlatfield();
        virtual void applyQuantumEfficiency();
        virtual void addPhotonNoise();
        virtual void applyFullWellSaturation();
        virtual void applyChargeTranferSmearing();
        virtual void addReadoutNoise();
        virtual void applyGain();
        virtual void addBias();
        virtual void applyDigitalSaturation();

        unsigned int Nrows;                      // Number of rows
        unsigned int Ncols;                      // Number of columns, readout is in column direction.
        unsigned int Nsubpixels;                 // Each pixel consists of Nsubpixels*Nsubpixels subpixels
        double readoutTime;                      // [s]
        double chargeTransferTime;               // [s]
        unsigned int biasLevel;                  // Electronic offset [ADU]
        unsigned int readoutNoise;               // [electrons]
        unsigned int gain;                       // [electrons/ADU]
        unsigned long digitalSaturationLimit;    // [ADU/pix]
        unsigned long fullWellCapacity;          // [electrons]
        double quantumEfficiency;                // Detected photon to converted electron ratio, in [0,1]
        double pixelSize;                        // [micrometer]
        double pixelScale;                       // [arcsec/pix]
        double originOffsetX;                    // x-coord of pixel (0,0) in focal plane reference system [mm]
        double originOffsetY;                    // y-coord of pixel (0,0) in focal plane reference system [mm] 

        double **pixelMap;                       // Par of a CCD image. Nrows x Ncols
 
        Telescope telescope;

        double internalTime;


    private:

};



#endif
