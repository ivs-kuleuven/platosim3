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


class Clock;

class Observer
{
    public:
        virtual ~Observer();
        virtual void update(ParaJitter* newJitter, double currentJitterStep) = 0;

    protected:
        Observer();
};

class Clock
{
	public:

		Clock();
		~Clock();

		void startSimulation();
		void updateJitter();

		void attach(Observer*);
		void detach(Observer*);
		void notify();

	protected:

        bool endSimulation;

    private:

    	std::vector<Observer*> observers;

}
