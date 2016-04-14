
#include "Logger.h"



// Logger::Logger()
// 
// PURPOSE: default constructor
//      
// INPUT: None
//
// OUTPUT: None
//

Logger::Logger()
{
    enabledLogLevels = DEBUG | INFO | WARNING | ERROR;
}








// Logger::addOutputStream()
// 
// PURPOSE: add an existing output stream to which log messages will be sent.
//		
// INPUT: logLevel: one of the log levels defined in Logger.h
//
// OUTPUT: None
//

void Logger::addOutputStream(ostream &outputStream, LogLevel logLevel)
{
	outputStreams.push_back(&outputStream);
	outputStreamLogLevel.push_back(logLevel);
}







// Logger:emit()
//
// PURPOSE: send a message to all logging output streams (e.g. cerr, cout, file stream, ...)
//          for which the given logLevel was set when the streams were added to the logger.
//
// INPUT: message: a string containing the log message
//        logLevel: one of the log levels defined in Logger.h
//
// OUTPUT: given 
//

void Logger::emit(string message, LogLevel logLevel)
{
    // Get the current time [clock ticks]

    auto currentTime = chrono::system_clock::to_time_t(chrono::system_clock::now());
    
    // Convert this into the date and time in the local time zone

    struct tm *localDateAndTime = localtime(&currentTime);
    
    // Format the result in a time stamp containing date and time

    char timeStamp[80];
    strftime(timeStamp, 80, "%F %T", localDateAndTime);

	if (outputStreams.size() != 0)
	{
    	for (int n = 0; n < outputStreams.size(); ++n)
    	{
    		if (outputStreamLogLevel[n] & logLevel & enabledLogLevels)
    		{
                *outputStreams[n] << timeStamp;
                *outputStreams[n] << " " << setw(7) << left << logLevelName[logLevel];
                *outputStreams[n] << " " << message << endl;
                (*outputStreams[n]).flush();
    		}
		}
	}
}





// Logger:debug()
//
// PURPOSE: a convenient alias for emit(message, DEBUG)
// 
// INPUT: message: see emit()
//
// OUTPUT: None

void Logger::debug(string message)
{
	emit(message, DEBUG);
}






// Logger::info()
//
// PURPOSE: a convenient alias for emit(message, INFO)
// 
// INPUT: message: see emit()
//
// OUTPUT: None

void Logger::info(string message)
{
	emit(message, INFO);
}





// Logger::warning()
//
// PURPOSE: a convenient alias for emit(message, WARNING)
// 
// INPUT: message: see emit()
//
// OUTPUT: None

void Logger::warning(string message)
{
	emit(message, WARNING);
}






// Logger::error()
//
// PURPOSE: a convenient alias for emit(message, ERROR)
// 
// INPUT: message: see emit()
//
// OUTPUT: None

void Logger::error(string message)
{
	emit(message, ERROR);
}







// Logger::enableLogLevel()
//
// PURPOSE: Make sure that sending log messages with the given logLevel are not ignored.
//          Can be undone with disableLogLevel(). 
// 
// INPUT: logLevel: one of the log levels specified in Logger.h
//
// OUTPUT: None

void Logger::enableLogLevel(LogLevel logLevel)
{
    enabledLogLevels = enabledLogLevels | logLevel;
}








// Logger::disableLogLevel()
//
// PURPOSE: log messages with the given logLevel are ignored.
//          For example, if disableLogLevel(DEBUG) was issued, sending log messages
//          with debug() are ignored. Can be undone with enableLogLevel()
//       
// 
// INPUT: logLevel: one of the log levels specified in Logger.h
//
// OUTPUT: None

void Logger::disableLogLevel(LogLevel logLevel)
{
    enabledLogLevels = enabledLogLevels & (~logLevel);
}

