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

	/**
	 * Start the clock.
	 * Returns the time at which the instance was started.
	 */
	my_steady_clock::time_point start() noexcept;

	/**
	 * Stops the clock permanently for the instance of the Timer.
	 * Returns the time at which the instance was stopped.
	 */
	my_steady_clock::time_point stop() noexcept;

	/**
	 * Returns a chrono::duration since the Timer was started.
	 * If stop() was called, the duration is between the start() and stop()
	 * calls is returned.
	 */
	std::chrono::duration<uint64_t, std::milli> elapsed() const noexcept;

	/**
	 * Returns the number of milliseconds since the Timer was started.
	 * If stop() was called, the duration is between the start() and stop()
	 * calls is returned.
	 */
	uint64_t elapsedMilliseconds() const noexcept;

	/**
	 * Returns the number of microseconds since the Timer was started.
	 * If stop() was called, the duration is between the start() and stop()
	 * calls is returned.
	 */
	uint64_t elapsedMicroseconds() const noexcept;

	/**
	 * Returns the number of nanoseconds since the Timer was started.
	 * If stop() was called, the duration is between the start() and stop()
	 * calls is returned.
	 */
	uint64_t elapsedNanoseconds() const noexcept;

	/**
	 * Returns the time at which the instance was started.
	 */
	my_steady_clock::time_point startTime() const noexcept;

	/**
	 * Returns the time at which the instance was stopped.
	 */
	my_steady_clock::time_point stopTime() const noexcept;

	/**
	 * Returns a human-readable representation including the elapsed time and unit.
	 */
	std::string elapsedTag() const;

protected:
	bool running{false};                   //!< true if timer has been started and not stopped after that
	my_steady_clock::time_point started;   //!< time at which timer has been started
	my_steady_clock::time_point stopped;   //!< time at which timer has been stopped

	/// If running returns now, otherwise the stop time
	my_steady_clock::time_point stopTimeOrNow() const noexcept;
};

/**
 * A timer for running time measurements.
 * Same as Timer but automatically starts on construction.
 */
class StartedTimer : public Timer {
public:
	StartedTimer() : Timer() {
		start();
	}
};

/**
 * A timer for running time measurements within a scope.
 *
 * Same as Timer but automatically starts on construction and report on destruction.
 *
 * @code
 * {
 *    Aux::ScopedTimer someName("Algorithm B");
 *    // some expensive operations
 *
 *
 * } // <- Report time since creation of the timer object to std::cout
 * @endcode
 *
 * @warning Similarly as std::unique_lock the timer needs a name ("someName" in the example)
 * to exist until the end of the scope. If the time measured seems too short, make you created
 * a named instances.
 */
class ScopedTimer : public StartedTimer {
public:
	explicit ScopedTimer(const std::string& label = "") :
		StartedTimer(),
		label(label)
	{}

	~ScopedTimer();

private:
	std::string label;
};

} /* namespace Aux */

#endif /* TIMER_H_ */

