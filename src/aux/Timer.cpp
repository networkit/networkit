/*
 * Timer.cpp
 *
 *  Created on: 14.01.2013
 *      Author: cls
 */

#include "Timer.h"

namespace Aux {

Timer::Timer() {
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

std::chrono::duration<int64_t, std::milli> Timer::elapsed() {
	auto now = std::chrono::steady_clock::now();
	// FIXME: return now - this->started;
}

std::chrono::steady_clock::time_point Timer::startTime() {
}

std::chrono::steady_clock::time_point Timer::stopTime() {
}

} /* namespace Aux */
