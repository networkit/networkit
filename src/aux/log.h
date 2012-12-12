/*
 * log.h
 *
 *  Created on: 28.11.2012
 *      Author: cls
 */

#ifndef LOG_H_
#define LOG_H_

#include "log4cxx/logger.h"

// prepends the function name
#define LOCATION "in " << __PRETTY_FUNCTION__ << ": "
#define LOGGER log4cxx::Logger::getRootLogger()


// short macros for logging statements
#define FATAL(X) LOG4CXX_FATAL(LOGGER, LOCATION << X)
#define ERROR(X) LOG4CXX_ERROR(LOGGER, LOCATION << X)
#define WARN(X) LOG4CXX_WARN(LOGGER, LOCATION << X)
#define INFO(X) LOG4CXX_INFO(LOGGER, LOCATION << X)
#define DEBUG(X) LOG4CXX_DEBUG(LOGGER, LOCATION << X);
#define TRACE(X) LOG4CXX_TRACE(LOGGER, LOCATION << X)


#endif /* LOG_H_ */
