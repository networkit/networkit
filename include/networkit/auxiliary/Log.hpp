#ifndef NETWORKIT_AUXILIARY_LOG_HPP_
#define NETWORKIT_AUXILIARY_LOG_HPP_

#include <sstream>
#include <string>

#include <networkit/auxiliary/StringBuilder.hpp>

#ifdef _MSC_VER
#define NETWORKT_PRETTY_FUNCTION __FUNCSIG__
#else
#define NETWORKT_PRETTY_FUNCTION __PRETTY_FUNCTION__
#endif // _MSC_VER

/// Logging without format string
#define LOG_AT(level, ...) \
    ::Aux::Log::log({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__}, level, __VA_ARGS__)
#define FATAL(...) LOG_AT(::Aux::Log::LogLevel::FATAL, __VA_ARGS__)
#define ERROR(...) LOG_AT(::Aux::Log::LogLevel::ERROR, __VA_ARGS__)
#define WARN(...) LOG_AT(::Aux::Log::LogLevel::WARN, __VA_ARGS__)
#define INFO(...) LOG_AT(::Aux::Log::LogLevel::INFO, __VA_ARGS__)

/// Logging with format string
#define LOG_ATF(level, ...) \
    ::Aux::Log::logF({__FILE__, NETWORKT_PRETTY_FUNCTION, __LINE__}, level, __VA_ARGS__)
#define FATALF(...) LOG_ATF(::Aux::Log::LogLevel::FATAL, __VA_ARGS__)
#define ERRORF(...) LOG_ATF(::Aux::Log::LogLevel::ERROR, __VA_ARGS__)
#define WARNF(...) LOG_ATF(::Aux::Log::LogLevel::WARN, __VA_ARGS__)
#define INFOF(...) LOG_ATF(::Aux::Log::LogLevel::INFO, __VA_ARGS__)

// DEBUG and TRACE are no-ops if NETWORKIT_RELEASE_LOGGING is defined.
#if defined(NETWORKIT_RELEASE_LOGGING)
#define NETWORKIT_LOG_DUMMY(...) \
    do { \
    } while (false)
#define DEBUG(...) NETWORKIT_LOG_DUMMY(__VA_ARGS__)
#define DEBUGF(...) NETWORKIT_LOG_DUMMY(__VA_ARGS__)
#define TRACE(...) NETWORKIT_LOG_DUMMY(__VA_ARGS__)
#define TRACEF(...) NETWORKIT_LOG_DUMMY(__VA_ARGS__)
#define TRACEPOINT NETWORKIT_LOG_DUMMY()
#else
#define DEBUG(...) LOG_AT(::Aux::Log::LogLevel::DEBUG, __VA_ARGS__)
#define TRACE(...) LOG_AT(::Aux::Log::LogLevel::TRACE, __VA_ARGS__)

#define DEBUGF(...) LOG_ATF(::Aux::Log::LogLevel::DEBUG, __VA_ARGS__)
#define TRACEF(...) LOG_ATF(::Aux::Log::LogLevel::TRACE, __VA_ARGS__)

#define TRACEPOINT LOG_AT(::Aux::Log::LogLevel::TRACE, "tracepoint")
#endif // defined(NETWORKIT_RELEASE_LOGGING)

namespace Aux {
namespace Log {

struct Location {
    const char *file;
    const char *function;
    const int line;
};

enum class LogLevel {
    TRACE,
    DEBUG,
    INFO,
    WARN,
    ERROR,
    FATAL,
    QUIET,         // Emits no log messages at all.
    trace = TRACE, // this + following added for backwards compatibility
    debug = DEBUG,
    info = INFO,
    warn = WARN,
    error = ERROR,
    fatal = FATAL,
    quiet = QUIET
};

/**
 * Accept loglevel as string and set.
 * @param logLevel as string
 */
void setLogLevel(const std::string &logLevel);

/**
 * @return current loglevel as string
 */
std::string getLogLevel();

namespace Impl {
void log(const Location &loc, LogLevel p, const std::string &msg);
} // namespace Impl

///! Returns true iff logging at the provided level is currently activated
bool isLogLevelEnabled(LogLevel p) noexcept;

template <typename... T>
void log(const Location &loc, LogLevel p, const T &...args) {
    if (!isLogLevelEnabled(p))
        return;

    std::stringstream stream;
    printToStream(stream, args...);
    Impl::log(loc, p, stream.str());
}

template <typename... T>
void logF(const Location &loc, LogLevel p, const std::string &format, const T &...args) {
    if (!isLogLevelEnabled(p))
        return;

    std::stringstream stream;
    printToStreamF(stream, format, args...);
    Impl::log(loc, p, stream.str());
}

} // namespace Log
} // namespace Aux

#endif // NETWORKIT_AUXILIARY_LOG_HPP_
