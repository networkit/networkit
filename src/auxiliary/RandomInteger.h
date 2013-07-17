/*
 * RandomInteger.h
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef RANDOMINTEGER_H_
#define RANDOMINTEGER_H_

#include <random>
#include <ctime>

namespace Aux {

class RandomInteger {

protected:

	std::random_device randomDevice;
	std::default_random_engine randomEngine;
	std::uniform_int_distribution<int64_t> distribution;

public:

	RandomInteger();

	virtual ~RandomInteger();

	/**
	 * TODO: Christian, please specify if the bounds are in- or exclusive!
	 */
	int64_t generate(int64_t lower, int64_t upper);

	/**
	 * Generate random integer in closed interval [lower, upper]. Faster generation,
	 * but without large strength in terms of randomness;
	 */
	int64_t generateFast(int64_t lower, int64_t upper);
};

} /* namespace Aux */
#endif /* RANDOMINTEGER_H_ */
