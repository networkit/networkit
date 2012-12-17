/*
 * log.h
 *
 *  Created on: 28.11.2012
 *      Author: cls
 */

#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#include <map>
#include <unordered_map>

#include "log4cxx/logger.h"


// prepends the function name
#define LOCATION "in " << __PRETTY_FUNCTION__ << ": "
#define LOGGER log4cxx::Logger::getRootLogger()


// short macros for logging statements
#define FATAL(X) LOG4CXX_FATAL(LOGGER, LOCATION << X)
#define ERROR(X) LOG4CXX_ERROR(LOGGER, LOCATION << X)
#define WARN(X) LOG4CXX_WARN(LOGGER, LOCATION << X)
#define INFO(X) LOG4CXX_INFO(LOGGER, LOCATION << X)
#define DEBUG(X) LOG4CXX_DEBUG(LOGGER, LOCATION << X)
#define TRACE(X) LOG4CXX_TRACE(LOGGER, LOCATION << X)




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





#endif /* LOG_H_ */
