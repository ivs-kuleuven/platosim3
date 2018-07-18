#ifndef CLOCK_H
#define CLOCK_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <unistd.h>
#include <iostream>
#include <mutex>
#include <condition_variable>

#include <vector>

#include <thread>

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

        virtual JitterGenerator* getJitterInstance() = 0;
        virtual DriftGenerator* getDriftInstance() = 0;

        virtual bool isClient(){return false;};

        virtual bool simulationEnd(){return true;};

    protected:
        Observer();
};

class Clock
{
	public:

		Clock(string inputFilename);
		~Clock();

		void startSimulation(std::condition_variable* cond_var, bool* notified, bool* newStep, std::mutex* m);
		//bool waitForNextStep();

		void attach(Observer*);
		void detach(Observer*);
		void notify();

		virtual void configure(ConfigurationParameters &configParam);

	protected:

        bool endSimulation;

    private:

     	double simulationTime;
     	double timeStep;

    	std::vector<Observer*> observers;

    	JitterGenerator *jitterInstance;
    	DriftGenerator *driftInstance;

    	double exposureTime; 
        double readoutTime; 
	    int numExposures;
};

#endif