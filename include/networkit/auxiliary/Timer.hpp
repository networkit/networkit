/*
 * Timer.h
 *
 *  Created on: 14.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <cstdint>
#include <string>

namespace Aux {

/**
 * A timer for running time measurements.
 */
class Timer {
public:
	#ifdef __MIC__
		using my_steady_clock = std::chrono::monotonic_clock;
	#else
		using my_steady_clock = std::chrono::steady_clock;
	#endif

	Timer() = default;
	virtual ~Timer() = default;

	/**
	 * Start the clock.
	 * Returns the time at which the instance was started.
	 */
	virtual my_steady_clock::time_point start() noexcept;

	/**
	 * Stops the clock permanently for the instance of the Timer.
	 * Returns the time at which the instance was stopped.
	 */
	virtual my_steady_clock::time_point stop() noexcept;

	/**
	 * The number of milliseconds since the current time that the Timer
	 * object was created.  If stop() was called, it is the number
	 * of seconds from the instance creation until stop() was called.
	 */
	virtual std::chrono::duration<uint64_t, std::milli> elapsed() const noexcept;

	/**
	 * The number of milliseconds since the current time that the Timer
	 * object was created. If stop() was called, it is the number
	 * of seconds from the instance creation until stop() was called.
	 */
	virtual uint64_t elapsedMilliseconds() const noexcept;

	/**
	 * The number of microseconds since the current time that the Timer
	 * object was created. If stop() was called, it is the number
	 * of seconds from the instance creation until stop() was called.
	 */
	virtual uint64_t elapsedMicroseconds() const noexcept;

	/**
	 * The number of nanoseconds since the current time that the Timer
	 * object was created. If stop() was called, it is the number
	 * of seconds from the instance creation until stop() was called.
	 */
	virtual uint64_t elapsedNanoseconds() const noexcept;

	/**
	 * Returns the time at which the instance was started.
	 */
	virtual my_steady_clock::time_point startTime() const noexcept;

	/**
	 * Returns the time at which the instance was stopped.
	 */
	virtual my_steady_clock::time_point stopTime() const noexcept;

	/**
	 * Returns a human-readable representation including the elapsed time and unit.
	 */
	virtual std::string elapsedTag() const;

protected:
	bool running{false};                   //!< true if timer has been started and not stopped after that
	my_steady_clock::time_point started;   //!< time at which timer has been started
	my_steady_clock::time_point stopped;   //!< time at which timer has been stopped

	/// If running returns now, otherwise the stop time
	my_steady_clock::time_point stopTimeOrNow() const noexcept;
};

} /* namespace Aux */
#endif /* TIMER_H_ */

