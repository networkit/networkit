/*
 * log.h
 *
 *  Created on: 28.11.2012
 *      Author: cls
 */

#ifndef LOG_H_
#define LOG_H_

#include "log4cxx/logger.h"

#define FATAL(X) LOG4CXX_FATAL(log4cxx::Logger::getRootLogger(), X)
#define ERROR(X) LOG4CXX_ERROR(log4cxx::Logger::getRootLogger(), X)
#define WARN(X) LOG4CXX_WARN(log4cxx::Logger::getRootLogger(), X)
#define INFO(X) LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), X)
#define DEBUG(X) LOG4CXX_DEBUG(log4cxx::Logger::getRootLogger(), X);
#define TRACE(X) LOG4CXX_TRACE(log4cxx::Logger::getRootLogger(), X)


#endif /* LOG_H_ */
