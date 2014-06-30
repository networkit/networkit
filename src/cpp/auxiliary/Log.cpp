#include <atomic>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <ios>

#include "Log.h"

namespace Aux { namespace Log {

void setLogLevel(std::string logLevel) {
	if (logLevel == "TRACE") {
		Settings::setLogLevel(LogLevel::trace);
	} else if (logLevel == "DEBUG") {
		Settings::setLogLevel(LogLevel::debug);
	} else if (logLevel == "INFO") {
		Settings::setLogLevel(LogLevel::info);
	} else if (logLevel == "WARN") {
		Settings::setLogLevel(LogLevel::warn);
	} else if (logLevel == "ERROR") {
		Settings::setLogLevel(LogLevel::error);		
	} else if (logLevel == "FATAL") {
		Settings::setLogLevel(LogLevel::fatal);
	} else {
		throw std::runtime_error("unknown loglevel");
	}
}

std::string getLogLevel() {
	LogLevel current = Settings::getLogLevel();
	if (current == LogLevel::trace) {
		return "TRACE";
	} else if (current == LogLevel::debug) {
		return "DEBUG";
	} else if (current == LogLevel::info) {
		return "INFO";
	} else if (current == LogLevel::warn) {
		return "WARN";
	} else if (current == LogLevel::error) {
		return "ERROR";
	} else if (current == LogLevel::fatal) {
		return "FATAL";
	} else {
		// this only exists to silence a warning:
		// TODO: consider replacing it with __builtin_unreachable();
		throw std::logic_error{"invalid loglevel. This should NEVER happen"};
	}
}

namespace Settings {

namespace {
bool printTime = false;
bool printLocation = false;
LogLevel loglevel = LogLevel::info;
std::ofstream logfile;

std::atomic_bool logfileIsOpen{false};
std::mutex logfileMutex;
}

LogLevel getLogLevel() {return loglevel;}
void setLogLevel(LogLevel p) {loglevel = p;}

void setPrintTime(bool b) {printTime = b;}
bool getPrintTime() {return printTime;}

void setPrintLocation(bool b) {printLocation = b;}
bool getPrintLocation() {return printLocation;}

void setLogfile(const std::string& filename) {
	std::lock_guard<std::mutex> guard{logfileMutex};
	if(logfile.is_open()) {
		logfile.close();
	}
	if(filename.empty()) {
		logfile.open(filename, std::ios_base::out | std::ios_base::app);
		logfileIsOpen = logfile.is_open();
	} else {
		logfileIsOpen = false;
	}
}
} // namespace Settings

void printLogLevel(std::ostream& stream, LogLevel p) {
	switch(p) {
		case LogLevel::fatal:
			stream << "[FATAL]"; break;
		case LogLevel::error:
			stream << "[ERROR]"; break;
		case LogLevel::warn:
			stream << "[WARN ]"; break;
		case LogLevel::info:
			stream << "[INFO ]"; break;
		case LogLevel::debug:
			stream << "[DEBUG]"; break;
		case LogLevel::trace:
			stream << "[TRACE]"; break;
	}
}

void printTime(std::ostream& stream,
		const std::chrono::time_point<std::chrono::system_clock>& timePoint) {
	stream << "[" << timePoint.time_since_epoch().count() << "]";
}

void printLocation(std::ostream& stream, const Location& loc) {
	stream << "[“" << loc.file << "”, " << loc.line << ": " << loc.function << "]";
}

std::tuple<std::string, std::string> getTerminalFormat(LogLevel p) {
	switch (p) {
		case LogLevel::fatal:
			return std::make_tuple("\033[1;31m", "\033[0m");
		case LogLevel::error:
			return std::make_tuple("\033[31m", "\033[0m");
		case LogLevel::warn:
			return std::make_tuple("\033[33m", "\033[0m");
		case LogLevel::info:
		case LogLevel::debug:
		case LogLevel::trace:
			return std::make_tuple("", "");
		default:
			// this only exists to silence a warning:
			// TODO: consider replacing it with __builtin_unreachable();
			throw std::logic_error{"invalid loglevel. This should NEVER happen"};
	}
}

static void logToTerminal(const Location& loc, LogLevel p,
		const std::chrono::time_point<std::chrono::system_clock>& timePoint,
		const std::string msg) {
	std::stringstream stream;
	
	if(Settings::getPrintTime()) {
		printTime(stream, timePoint);
	}
	
	std::string termFormatOpen, termFormatClose;
	std::tie(termFormatOpen, termFormatClose) = getTerminalFormat(p);
	
	stream << termFormatOpen;
	printLogLevel(stream, p);
	stream <<termFormatClose;
	
	if(Settings::getPrintLocation()) {
		printLocation(stream, loc);
	}
	
	stream << ": ";
	
	stream << termFormatOpen;
	stream << msg;
	stream << termFormatClose;
	
	stream.put('\n');
	
	if (p < LogLevel::warn) {
		static std::mutex cout_mutex;
		{
			std::lock_guard<std::mutex> guard{cout_mutex};
			std::cout << stream.str();
		}
	} else {
		static std::mutex cerr_mutex;
		{
			std::lock_guard<std::mutex> guard{cerr_mutex};
			std::cerr << stream.str();
		}
	}
}

static void logToFile(const Location& loc, LogLevel p,
		const std::chrono::time_point<std::chrono::system_clock>& timePoint,
		const std::string& msg) {
	if(!Settings::logfileIsOpen) {
		return;
	}
	std::stringstream stream;
	printTime(stream, timePoint);
	stream << ' ';
	printLogLevel(stream, p);
	
	if(Settings::getPrintLocation()) {
		stream << ' ';
		printLocation(stream, loc);
	}
	
	stream << ": " << msg << '\n';
	{
		std::lock_guard<std::mutex> guard{Settings::logfileMutex};
		if(!Settings::logfileIsOpen) {
			return;
		}
		Settings::logfile << stream.str() << std::flush;
	}
}

namespace Impl {

void log(const Location& loc, LogLevel p, const std::string msg) {
	auto time =std::chrono::system_clock::now();
	
	logToTerminal(loc, p, time, msg);
	logToFile(loc, p, time, msg);
}

} // namespace impl

}} // namespace Aux::Log
