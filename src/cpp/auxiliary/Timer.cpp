/*
 * Timer.cpp
 *
 *  Created on: 14.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "Timer.h"

namespace Aux {

Timer::Timer() : running(false) {
	// TODO Auto-generated constructor stub

}

Timer::~Timer() {
	// TODO Auto-generated destructor stub
}

std::chrono::steady_clock::time_point Timer::start() {
	this->started = std::chrono::steady_clock::now();
	return this->started;
}

std::chrono::steady_clock::time_point Timer::stop() {
	this->stopped = std::chrono::steady_clock::now();
	return this->stopped;
}

std::chrono::duration<uint64_t, std::milli> Timer::elapsed() {
	std::chrono::duration<uint64_t, std::milli> elapsed = std::chrono::duration_cast<std::chrono::duration<uint64_t, std::milli>>(this->stopped - this->started);
	return elapsed;
}

std::chrono::steady_clock::time_point Timer::startTime() {
	return this->started;
}

std::chrono::steady_clock::time_point Timer::stopTime() {
	return this->stopped;
}

uint64_t Timer::elapsedMilliseconds() {
	return this->elapsed().count();
}

std::string Timer::elapsedTag() {
	std::stringstream s;
	s << "(" << this->elapsed().count() << " ms) ";
	return s.str();
}


} /* namespace Aux */

