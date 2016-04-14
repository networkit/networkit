#ifndef LOG_H_
#define LOG_H_

#include <iostream>
#include <mutex>

#include "StringBuilder.h"

#ifdef NOLOGGING

#define FATAL(...) do{}while(false)
#define ERROR(...) do{}while(false)
#define WARN(...)  do{}while(false)
#define INFO(...)  do{}while(false)
#define DEBUG(...) do{}while(false)
#define TRACE(...) do{}while(false)

#define FATALF(...) do{}while(false)
#define ERRORF(...) do{}while(false)
#define WARNF(...)  do{}while(false)
#define INFOF(...)  do{}while(false)
#define DEBUGF(...) do{}while(false)
#define TRACEF(...) do{}while(false)

#define TRACEPOINT do{}while(false)

#else // NOLOGGING

#define LOG_LEVEL_FATAL 0
#define LOG_LEVEL_ERROR 1
#define LOG_LEVEL_WARN  2
#define LOG_LEVEL_INFO  3
#define LOG_LEVEL_DEBUG 4
#define LOG_LEVEL_TRACE 5


#if !defined(LOG_LEVEL)
#define LOG_LEVEL LOG_LEVEL_TRACE
#endif

#if LOG_LEVEL >= LOG_LEVEL_FATAL
#define FATAL(...) ::Aux::Log::log({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::fatal, __VA_ARGS__)
#define FATALF(...) ::Aux::Log::logF({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::fatal, __VA_ARGS__)
#else
#define FATAL(...) do{}while(false)
#define FATALF(...) do{}while(false)
#endif

#if LOG_LEVEL >= LOG_LEVEL_ERROR
#define ERROR(...) ::Aux::Log::log({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::error, __VA_ARGS__)
#define ERRORF(...) ::Aux::Log::logF({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::error, __VA_ARGS__)
#else
#define ERROR(...) do{}while(false)
#define ERRORF(...) do{}while(false)
#endif

#if LOG_LEVEL >= LOG_LEVEL_WARN
#define WARN(...)  ::Aux::Log::log({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::warn,  __VA_ARGS__)
#define WARNF(...)  ::Aux::Log::logF({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::warn,  __VA_ARGS__)
#else
#define WARN(...) do{}while(false)
#define WARNF(...) do{}while(false)
#endif

#if LOG_LEVEL >= LOG_LEVEL_INFO
#define INFO(...)  ::Aux::Log::log({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::info,  __VA_ARGS__)
#define INFOF(...)  ::Aux::Log::logF({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::info,  __VA_ARGS__)
#else
#define INFO(...) do{}while(false)
#define INFOF(...) do{}while(false)
#endif

#if LOG_LEVEL >= LOG_LEVEL_DEBUG
#define DEBUG(...) ::Aux::Log::log({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::debug, __VA_ARGS__)
#define DEBUGF(...) ::Aux::Log::logF({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::debug, __VA_ARGS__)
#else
#define DEBUG(...) do{}while(false)
#define DEBUGF(...) do{}while(false)
#endif

#if LOG_LEVEL >= LOG_LEVEL_TRACE
#define TRACE(...) ::Aux::Log::log({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::trace, __VA_ARGS__)
#define TRACEF(...) ::Aux::Log::logF({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::trace, __VA_ARGS__)
#define TRACEPOINT ::Aux::Log::log({__FILE__, __PRETTY_FUNCTION__, __LINE__},\
		::Aux::Log::LogLevel::trace, "tracepoint")
#else
#define TRACE(...) do{}while(false)
#define TRACEF(...) do{}while(false)
#define TRACEPOINT do{}while(false)
#endif

#endif // NOLOGGING

namespace Aux { namespace Log {

struct Location {
	const char* file;
	const char* function;
	const int line;
};

enum class LogLevel {
	trace,
	debug,
	info,
	warn,
	error,
	fatal
};

/**
 * Accept loglevel as string and set.
 * @param logLevel as string
 */
void setLogLevel(std::string logLevel);

/**
 * @return current loglevel as string
 */
std::string getLogLevel();

namespace Settings {

LogLevel getLogLevel();
void setLogLevel(LogLevel p = LogLevel::info);

void setPrintTime(bool b);
bool getPrintTime();

void setPrintLocation(bool b);
bool getPrintLocation();

void setLogfile(const std::string& filename);
}

namespace Impl {
void log(const Location& loc, LogLevel p, const std::string msg);
} //namespace impl

template<typename...T>
void log(const Location& loc, LogLevel p, const T&...args) {
	if(p >= Settings::getLogLevel()) {
		std::stringstream stream;
		printToStream(stream, args...);
		Impl::log(loc, p, stream.str());
	}
}

template<typename...T>
void logF(const Location& loc, LogLevel p, const std::string& format, const T&...args) {
	if(p >= Settings::getLogLevel()) {
		std::stringstream stream;
		printToStreamF(stream, format, args...);
		Impl::log(loc, p, stream.str());
	}
}


}} // namespace Aux::Log

#endif
