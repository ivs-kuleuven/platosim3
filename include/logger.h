
#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <ostream>
#include <iomanip>
#include <vector>
#include <map>
#include <chrono>
#include <ctime>


using namespace std;


// Define the the different log levels.

typedef unsigned short LogLevel;

static LogLevel DEBUG   = 1;
static LogLevel INFO    = 2;
static LogLevel WARNING = 4;
static LogLevel ERROR   = 8; 

static map<LogLevel, string> logLevelName{{DEBUG, "DEBUG"}, {INFO, "INFO"}, {WARNING, "WARNING"}, {ERROR, "ERROR"}};


// Define the Logger class

class Logger
{
	public:

		Logger();
		~Logger(){};

		void addOutputStream(ostream &newOutputStream, LogLevel logLevel);
		void debug(string message);
		void info(string message);
		void warning(string message);
		void error(string message);
		void enableLogLevel(LogLevel logLevel);
		void disableLogLevel(LogLevel logLevel);

	protected:

		void emit(string message, LogLevel logLevel);

	private:

		vector<ostream*> outputStreams;
		vector<LogLevel> outputStreamLogLevel;
		LogLevel enabledLogLevels;

};


// The Logger object is defined in main.cpp
// Every file including this header should know that the Logger
// is externally defined.

extern Logger Log;


#endif