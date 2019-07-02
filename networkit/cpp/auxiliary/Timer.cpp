/*
 * Timer.cpp
 *
 *  Created on: 14.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <networkit/auxiliary/Timer.hpp>

namespace Aux {

Timer::my_steady_clock::time_point Timer::start() {
	started = my_steady_clock::now();
	running = true;
	return started;
}

Timer::my_steady_clock::time_point Timer::stop() {
	stopped = my_steady_clock::now();
	running = false;
	return stopped;
}

std::chrono::duration<uint64_t, std::milli> Timer::elapsed() const {
	if (running) {
		return std::chrono::duration_cast<std::chrono::duration<uint64_t, std::milli>>(std::chrono::steady_clock::now() - this->started);
	}
	std::chrono::duration<uint64_t, std::milli> elapsed = std::chrono::duration_cast<std::chrono::duration<uint64_t, std::milli>>(this->stopped - this->started);
	return elapsed;
}

Timer::my_steady_clock::time_point Timer::startTime() {
	return started;
}

Timer::my_steady_clock::time_point Timer::stopTime() {
	return stopped;
}

uint64_t Timer::elapsedMilliseconds() {
	return elapsed().count();
}

uint64_t Timer::elapsedMicroseconds() {
	return std::chrono::duration_cast<std::chrono::duration<uint64_t, std::micro>>(this->stopped - this->started).count();
}

uint64_t Timer::elapsedNanoseconds() {
	return std::chrono::duration_cast<std::chrono::duration<uint64_t, std::nano>>(this->stopped - this->started).count();
}

std::string Timer::elapsedTag() {
	std::stringstream s;
	s << "(" << elapsedMilliseconds() << " ms) ";
	return s.str();
}


} /* namespace Aux */
