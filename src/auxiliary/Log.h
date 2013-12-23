/*
 * log.h
 *
 *  Created on: 28.11.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <sstream>
#include "Loglevel.h"

namespace Aux {

// prepends the function name
#define LOCATION "in " << __PRETTY_FUNCTION__ << ": "


// short macros for logging statements
#if defined NOLOGGING
	// define logging macros to nothing
	#define FATAL(X)
	#define ERROR(X)
	#define WARN(X)
	#define INFO(X)
	#define DEBUG(X)
	#define TRACE(X)
#else
	#define LOG_TEST_LEVEL(v)   v >= Aux::LogLevel::getLevel()
	#define FATAL(X) if(LOG_TEST_LEVEL(5)) std::cout << "[FATAL] " << LOCATION << X << std::endl;
	#define ERROR(X) if(LOG_TEST_LEVEL(4)) std::cout << "[ERROR] " << LOCATION << X << std::endl;
	#define WARN(X) if(LOG_TEST_LEVEL(3)) std::cout << "[WARN] " << LOCATION << X << std::endl;
	#define INFO(X) if(LOG_TEST_LEVEL(2)) std::cout << "[INFO] " << LOCATION << X << std::endl;
	#define DEBUG(X) if(LOG_TEST_LEVEL(1)) std::cout << "[DEBUG] " << LOCATION << X << std::endl;
	#define TRACE(X) if(LOG_TEST_LEVEL(0)) std::cout << "[TRACE] " << LOCATION << X << std::endl;
/**
 * Set the current loglevel.
 * @param loglevel 
 */
inline void setLoglevel(const std::string& loglevel) {
	if (loglevel == "TRACE") {
		Aux::LogLevel::setLevel(0);
	} else if (loglevel == "DEBUG") {
		Aux::LogLevel::setLevel(1);
	} else if (loglevel == "INFO") {
		Aux::LogLevel::setLevel(2);
	} else if (loglevel == "WARN") {
		Aux::LogLevel::setLevel(3);
	} else if (loglevel == "ERROR") {
		Aux::LogLevel::setLevel(4);
	} else {
		ERROR("unknown loglevel: " << loglevel);
		exit(1);
	}
}


/**
 * Call this first to configure logging output.
 */
inline void configureLogging(const std::string& loglevel = "ERROR") {
	// configure logging
	//log4cxx::BasicConfigurator::configure();
	setLoglevel(loglevel);
}


/**
 * Get the current log level.
 */
inline std::string currentLogLevel() {
	std::string s;
	switch (Aux::LogLevel::getLevel()) {
		case 5: s = "FATAL"; break;
		case 4: s = "ERROR"; break;
		case 3: s = "WARN"; break;
		case 2: s = "INFO"; break;
		case 1: s = "DEBUG"; break;
		case 0: s = "TRACE"; break;
		default: s = "Loglevel is messed up, call setLoglevel([ERROR|WARN|INFO|DEBUG|TRACE])"; break;
	}
	return s;
}
#endif


#define PRINTMAP(M) std::cout << M << std::endl;
/**
 * Print a map for debugging.
 */
template <typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::map<K, V>& m)
{
    os << "{ ";
    for (auto& kv : m)
    {
        os << kv.first << ": " << kv.second << ", ";
    }
    return os << " }";
}

/**
 * Print an unordered_map for debugging.
 */
template <typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<K, V>& m)
{
    os << "{ ";
    for (typename std::unordered_map<K, V>::const_iterator i = m.begin(); i != m.end(); ++i)
    {
        if (i != m.begin()) os << ", ";
        os << i->first << ": " << i->second;
    }
    return os << " }";
}


} // end namespace



#endif /* LOG_H_ */
