#ifndef DETECTORFACTORY_H
#define DETECTORFACTORY_H

#include <string>
#include <cmath>
#include <random>
#include <functional>
#include <valarray>

#include "Detector.h"
#include "DetectorWithSymmetricalMappedPSF.h"
#include "DetectorWithAsymmetricalMappedPSF.h"
#include "DetectorWithAnalyticGaussianPSF.h"
#include "DetectorWithAnalyticNonGaussianPSF.h"

#include "ClosedLoopDetectorClasses.h"

using namespace std;

// using the abstract factory pattern to be able to create detector instances depending on whether a closed loop
// test is running or not

class AbstractDetectorFactory
{
    public:

        virtual Detector* createDetectorWithSymmetricalMappedPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure) = 0;

        virtual Detector* createDetectorWithAsymmetricalMappedPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure) = 0;

        virtual Detector* createDetectorWithAnalyticGaussianPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure) = 0;

        virtual Detector* createDetectorWithAnalyticNonGaussianPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure) = 0;
};

// this is the factory that creates the usual detector variants (Mapped, AnalyticNonGaussian, AnalyticGaussian) 

class DetectorFactory : public AbstractDetectorFactory
{
    public:

        Detector* createDetectorWithSymmetricalMappedPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
        {
            return new DetectorWithSymmetricalMappedPSF(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure);
        };

        Detector* createDetectorWithAsymmetricalMappedPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
        {
            return new DetectorWithAsymmetricalMappedPSF(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure);
        };

        Detector* createDetectorWithAnalyticGaussianPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
        {
            return new DetectorWithAnalyticGaussianPSF(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure);
        };

        Detector* createDetectorWithAnalyticNonGaussianPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
        {
            return new DetectorWithAnalyticNonGaussianPSF(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure);
        };
};


// this is the factory that creates the closed loop test detector variants (Mapped, AnalyticNonGaussian, AnalyticGaussian) 

class ClosedLoopDetectorFactory : public AbstractDetectorFactory
{
    public:

        Detector* createDetectorWithSymmetricalMappedPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
        {
            return new ClosedLoopDetectorWithSymmetricalMappedPSF(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure);
        };

        Detector* createDetectorWithAsymmetricalMappedPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
        {
            return new ClosedLoopDetectorWithAsymmetricalMappedPSF(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure);
        };

        Detector* createDetectorWithAnalyticGaussianPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
        {
            return new ClosedLoopDetectorWithAnalyticGaussianPSF(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure);
        };

        Detector* createDetectorWithAnalyticNonGaussianPsfInstance(ConfigurationParameters &configParam, HDF5File &hdf5file, Camera &camera, TemperatureGenerator &feeTemperatureGenerator, TemperatureGenerator &detectorTemperatureGenerator, double readoutTimeBeforeNextExposure, double readoutTimeDuringNextExposure)
        {
            return new ClosedLoopDetectorWithAnalyticNonGaussianPSF(configParam, hdf5file, camera, feeTemperatureGenerator, detectorTemperatureGenerator, readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure);
        };
};

#endif
