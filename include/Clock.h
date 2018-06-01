#ifndef CLOCK_H
#define CLOCK_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <unistd.h>
#include <iostream>

#include <vector>
#include <ctime>

#include <chrono>
#include <thread>

#include <time.h>

#include <omp.h>

#include "ConfigurationParameters.h"

#include "JitterGenerator.h"
#include "NoJitter.h"
#include "JitterFromFile.h"
#include "JitterFromRedNoise.h"
#include "JitterFromNetwork.h"
#include "DriftGenerator.h"
#include "NoDrift.h"
#include "ThermoElasticDriftFromFile.h"
#include "ThermoElasticDriftFromRedNoise.h"

class Clock;

class Observer
{
    public:
        virtual ~Observer();
        virtual void update(double simulationTime) = 0;

    protected:
        Observer();
};

class Clock
{
	public:

		Clock(string inputFilename);
		~Clock();

		void startSimulation();
		bool waitForNextStep();

		void attach(Observer*);
		void detach(Observer*);
		void notify();

		virtual void configure(ConfigurationParameters &configParam);

	protected:

        bool endSimulation;

    private:

     	double simulationTime;

    	std::vector<Observer*> observers;

    	JitterGenerator *jitterInstance;
    	DriftGenerator *driftInstance;

    	double exposureTime; 
	    int beginExposureNr;
	    int numExposures;
	    double readoutTime; 
};

#endif