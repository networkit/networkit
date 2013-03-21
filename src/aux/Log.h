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
#include <log4cxx/logger.h>
#endif

namespace Aux {

// prepends the function name
#define LOCATION "in " << __PRETTY_FUNCTION__ << ": "
#define LOGGER log4cxx::Logger::getRootLogger()


// short macros for logging statements
#ifdef NOLOG4CXX
	// define logging macros to nothing
	#define FATAL(X)
	#define ERROR(X)
	#define WARN(X)
	#define INFO(X)
	#define DEBUG(X)
	#define TRACE(X)
#else
	#define FATAL(X) LOG4CXX_FATAL(LOGGER, LOCATION << X)
	#define ERROR(X) LOG4CXX_ERROR(LOGGER, LOCATION << X)
	#define WARN(X) LOG4CXX_WARN(LOGGER, LOCATION << X)
	#define INFO(X) LOG4CXX_INFO(LOGGER, LOCATION << X)
	#define DEBUG(X) LOG4CXX_DEBUG(LOGGER, LOCATION << X)
	#define TRACE(X) LOG4CXX_TRACE(LOGGER, LOCATION << X)
#endif


/**
 * String representation of std::vector<T>
 * 		[x, y, ..., z]
 */
template <typename T> std::string vectorToString(const std::vector<T>& vec) {
	std::ostringstream ss;
	ss << '[';
	bool first = true;
	for (T element : vec) {
	    if (!first) {
	        ss << ", ";
	    }
	    ss << element;
	    first = false;
	}
	ss << ']';
	return ss.str();
}

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
