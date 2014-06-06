/*
 * ProgressMeter.h
 *
 *  Created on: 02.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef PROGRESSMETER_H_
#define PROGRESSMETER_H_

#include <iostream>

namespace Aux {

/**
 * A rudimentary console progress bar.
 */
class ProgressMeter {

protected:
	int64_t n;	//!< number of elements
	int64_t i;	//!< interval

public:

	ProgressMeter(int64_t n, int64_t i): n{n}, i{i} {}

	virtual ~ProgressMeter() = default;

	/**
	 * Send a signal to the progress meter - constructor parameters decide when output is produced.
	 */
	inline void signal(int64_t v) {
		if ((v % this->i) == 0) {
			std::cout << "." << std::flush;
		}
	}

	 template<typename Output> inline void signal(int64_t v, Output out) {
			if ((v % this->i) == 0) {
				std::string o = out();
				std::cout << o << std::flush;
			}
	}

	/**
	 * Call this when process finished.
	 */
	inline void end() {
		std::cout << std::endl;
	}


};

} /* namespace Aux */
#endif /* PROGRESSMETER_H_ */
