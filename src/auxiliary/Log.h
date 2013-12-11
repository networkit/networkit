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

#ifndef NOLOGGING
#ifndef NOLOG4CXX
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#endif
#endif

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
#elif defined SIMPLE
	#define LOG_TEST_CONDITION        output
	#define LOG_TEST_LEVEL(v)   (v) <= DebugCake::DebugLevel
	//TODO: rename semantically, D=> LOGSTREAM, struct => class
	//TODO: move to own class with header and source.
	template<typename C,typename T = std::char_traits<C> >
	struct DInfo {
	    typedef std::basic_ostream<C,T>& (*StrFunc)(std::basic_ostream<C,T>&);
	    static std::basic_ostream<C,T>& getDefault();
	};

	struct DebugCake {
	    static int DebugLevel;
	};
	template<typename C,typename T = std::char_traits<C> >
	struct D {
	    D(int level)
	        :stream(DInfo<C,T>::getDefault())
	        ,output( LOG_TEST_LEVEL(level) )
	    {}
	    D(int level, std::basic_ostream<C,T>& s)
	        :stream(s)
	        ,output( LOG_TEST_LEVEL(level) )
	    {}

	    template<typename V>
	    D& operator<<(V const& value) {
	        if (LOG_TEST_CONDITION) {
	            stream << value;
	        }
	        return *this;
	    }

	    D& operator<<(typename DInfo<C,T>::StrFunc func) {
	        if (LOG_TEST_CONDITION) {
	            func(stream);
        	}
        	return *this;
	    }
	    private:
	    std::basic_ostream<C,T>&  stream;
	    bool                      output;

	};

	template<>
	std::ostream&  DInfo<char,std::char_traits<char>       >::getDefault()
	{return std::cout; }

	template<>
	std::wostream& DInfo<wchar_t,std::char_traits<wchar_t> >::getDefault()
	{return std::wcout; }

	typedef D<char>    Debug;
	typedef D<char>    SIMPLELOG;
	typedef D<wchar_t> WDebug;
	int DebugCake::DebugLevel = 3;
	#define FATAL(X) Aux::D<char>(5) << "[FATAL] " << LOCATION << X << std::endl;
	#define ERROR(X) Aux::D<char>(4) << "[ERROR] " << LOCATION << X << std::endl;
	#define WARN(X) Aux::D<char>(3) << "[WARN] " << LOCATION << X << std::endl;
	#define INFO(X) Aux::D<char>(2) << "[INFO] " << LOCATION << X << std::endl;
	#define DEBUG(X) Aux::D<char>(1) << "[DEBUG] " << LOCATION << X << std::endl;
	#define TRACE(X) Aux::D<char>(0) << "[TRACE] " << LOCATION << X << std::endl;
/**
 * Set the current loglevel.
 * @param loglevel 
 */
inline void setLoglevel(const std::string& loglevel) {
	if (loglevel == "TRACE") {
		DebugCake::DebugLevel = 0;
	} else if (loglevel == "DEBUG") {
		DebugCake::DebugLevel = 1;
	} else if (loglevel == "INFO") {
		DebugCake::DebugLevel = 2;
	} else if (loglevel == "WARN") {
		DebugCake::DebugLevel = 3;
	} else if (loglevel == "ERROR") {
		DebugCake::DebugLevel = 4;
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
	switch (DebugCake::DebugLevel) {
		case 5: s = "FATAL"; break;
		case 4: s = "ERROR"; break;
		case 3: s = "WARNING"; break;
		case 2: s = "INFO"; break;
		case 1: s = "DEBUG"; break;
		case 0: s = "TRACE"; break;
		default: s = "Loglevel is messed up, call setLoglevel([ERROR|WARNING|INFO|DEBUG|TRACE])"; break;
	}
	return s;
}
#elif defined NOLOG4CXX
	// define logging macros to use std::cout
	#define FATAL(X) std::cout << "[FATAL] " << LOCATION << X << std::endl;
	#define ERROR(X) std::cout << "[ERROR] " << LOCATION << X << std::endl;
	#define WARN(X) std::cout << "[WARN] " << LOCATION << X << std::endl;
	#define INFO(X) std::cout << "[INFO] " << LOCATION << X << std::endl;
	#define DEBUG(X) std::cout << "[DEBUG] " << LOCATION << X << std::endl;
	#define TRACE(X) std::cout << "[TRACE] " << LOCATION << X << std::endl;
#else
	#define LOGGER log4cxx::Logger::getRootLogger()
	#define FATAL(X) LOG4CXX_FATAL(LOGGER, LOCATION << X)
	#define ERROR(X) LOG4CXX_ERROR(LOGGER, LOCATION << X)
	#define WARN(X) LOG4CXX_WARN(LOGGER, LOCATION << X)
	#define INFO(X) LOG4CXX_INFO(LOGGER, LOCATION << X)
	#define DEBUG(X) LOG4CXX_DEBUG(LOGGER, LOCATION << X)
	#define TRACE(X) LOG4CXX_TRACE(LOGGER, LOCATION << X)
#endif

#ifndef SIMPLE
#ifndef NOLOGGING
#ifndef NOLOG4CXX


/**
 * Set the current loglevel.
 * @param loglevel 
 */
inline void setLoglevel(const std::string& loglevel) {
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
 * Call this first to configure logging output.
 */
inline void configureLogging(const std::string& loglevel = "ERROR") {
	// configure logging
	log4cxx::BasicConfigurator::configure();
	setLoglevel(loglevel);
}


/**
 * Get the current log level.
 */
inline std::string currentLogLevel() {
	std::ostringstream s;
	s << log4cxx::Logger::getRootLogger()->getLevel()->toString();
	return s.str();
}
#endif //NOLOG4CXX
#endif //NOLOGGING
#endif //SIMPLE


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
