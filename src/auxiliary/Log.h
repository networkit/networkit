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

#ifndef NOLOG4CXX
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#endif

namespace Aux {

// prepends the function name
#define LOCATION "in " << __PRETTY_FUNCTION__ << ": "
#define LOGGER log4cxx::Logger::getRootLogger()


// short macros for logging statements
#if defined NOLOGGING
	// define logging macros to nothing
	#define FATAL(X)
	#define ERROR(X)
	#define WARN(X)
	#define INFO(X)
	#define DEBUG(X)
	#define TRACE(X)
#elif defined NOLOG4CXX
	// define logging macros to use std::cout
	#define FATAL(X) std::cout << "[FATAL] " << LOCATION << X << std::endl;
	#define ERROR(X) std::cout << "[ERROR] " << LOCATION << X << std::endl;
	#define WARN(X) std::cout << "[WARN] " << LOCATION << X << std::endl;
	#define INFO(X) std::cout << "[INFO] " << LOCATION << X << std::endl;
	#define DEBUG(X) std::cout << "[DEBUG] " << LOCATION << X << std::endl;
	#define TRACE(X) std::cout << "[TRACE] " << LOCATION << X << std::endl;
#else
	#define FATAL(X) LOG4CXX_FATAL(LOGGER, LOCATION << X)
	#define ERROR(X) LOG4CXX_ERROR(LOGGER, LOCATION << X)
	#define WARN(X) LOG4CXX_WARN(LOGGER, LOCATION << X)
	#define INFO(X) LOG4CXX_INFO(LOGGER, LOCATION << X)
	#define DEBUG(X) LOG4CXX_DEBUG(LOGGER, LOCATION << X)
	#define TRACE(X) LOG4CXX_TRACE(LOGGER, LOCATION << X)
#endif


#ifndef NOLOGGING
#ifndef NOLOG4CXX
/**
 * Call this first to configure logging output.
 */
inline void configureLogging(const std::string& loglevel = "ERROR") {
	// configure logging
	log4cxx::BasicConfigurator::configure();
	if (loglevel == "TRACE") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getTrace());
	} else if (loglevel == "DEBUG") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getDebug());
	} else if (loglevel == "INFO") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getInfo());
	} else if (loglevel == "WARN") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getWarn());
	} else if (loglevel == "ERROR") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getError());
	} else {
		ERROR("unknown loglevel: " << loglevel);
		exit(1);
	}
}

/**
 * Get the current log level.
 */
inline std::string currentLogLevel() {
	std::ostringstream s;
	s << log4cxx::Logger::getRootLogger()->getLevel()->toString();
	return s.str();
}
#endif
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
