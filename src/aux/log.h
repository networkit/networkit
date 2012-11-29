/*
 * log.h
 *
 *  Created on: 28.11.2012
 *      Author: cls
 */

#ifndef LOG_H_
#define LOG_H_


#include "log4cxx/logger.h"

// TODO: define concise logging macros here

#define DEBUG(x) LOG4CXX_DEBUG(log4cxx::Logger::getRootLogger(), x);
// FIXME: why doesn't this macro work?

#endif /* LOG_H_ */
