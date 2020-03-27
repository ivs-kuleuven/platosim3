#ifndef CLOSEDLOOPDETECTORCLASSES_H
#define CLOSEDLOOPDETECTORCLASSES_H

#include "Detector.h"
#include "DetectorWithSymmetricalMappedPSF.h"
#include "DetectorWithAsymmetricalMappedPSF.h"
#include "DetectorWithAnalyticGaussianPSF.h"
#include "DetectorWithAnalyticNonGaussianPSF.h"

#include "ClosedLoopUtility.h"

// variant of the detectorWithMappedPsf class

class ClosedLoopDetectorWithSymmetricalMappedPSF: public DetectorWithSymmetricalMappedPSF, public ClosedLoopUtility
{

    public:

        ClosedLoopDetectorWithSymmetricalMappedPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)       
        : DetectorWithSymmetricalMappedPSF{configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure}
        , ClosedLoopUtility{configParam}
        {}

        ~ClosedLoopDetectorWithSymmetricalMappedPSF() {};

    protected:

        virtual double takeExposure(int exposureNr, double startTime, double exposureTime) override;

        virtual void writePixelMapsToHDF5(int exposureNr) override;

        void setNewWindowPosition(std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple);

    private:


};


// variant of the detectorWithMappedPsf class

class ClosedLoopDetectorWithAsymmetricalMappedPSF: public DetectorWithAsymmetricalMappedPSF, public ClosedLoopUtility
{

    public:

        ClosedLoopDetectorWithAsymmetricalMappedPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)       
        : DetectorWithAsymmetricalMappedPSF{configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure}
        , ClosedLoopUtility{configParam}
        {}

        ~ClosedLoopDetectorWithAsymmetricalMappedPSF() {};

    protected:

        virtual double takeExposure(int exposureNr, double startTime, double exposureTime) override;

        virtual void writePixelMapsToHDF5(int exposureNr) override;

        void setNewWindowPosition(std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple);

    private:


};

// variant of the detectorWithAnalyticNonGaussian class

class ClosedLoopDetectorWithAnalyticNonGaussianPSF: public DetectorWithAnalyticNonGaussianPSF, public ClosedLoopUtility
{

    public:

        ClosedLoopDetectorWithAnalyticNonGaussianPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
        : DetectorWithAnalyticNonGaussianPSF{configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure}
        , ClosedLoopUtility{configParam}
        {}     

        ~ClosedLoopDetectorWithAnalyticNonGaussianPSF() {};

    protected:

        virtual double takeExposure(int exposureNr, double startTime, double exposureTime) override;

        virtual void writePixelMapsToHDF5(int exposureNr) override;

        void setNewWindowPosition(std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple);

    private:


};

// variant of the detectorWithAnalyticGaussian class

class ClosedLoopDetectorWithAnalyticGaussianPSF: public DetectorWithAnalyticGaussianPSF, public ClosedLoopUtility
{

    public:

        ClosedLoopDetectorWithAnalyticGaussianPSF(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)       
        : DetectorWithAnalyticGaussianPSF{configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure}
        , ClosedLoopUtility{configParam}
        {}

        ~ClosedLoopDetectorWithAnalyticGaussianPSF() {};

    protected:

        virtual double takeExposure(int exposureNr, double startTime, double exposureTime) override;

        virtual void writePixelMapsToHDF5(int exposureNr) override;

        void setNewWindowPosition(std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple);

    private:


};

#endif