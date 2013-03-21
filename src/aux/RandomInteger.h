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
	std::uniform_int_distribution<> distribution;

public:

	RandomInteger(int64_t lower, int64_t upper);

	virtual ~RandomInteger();

	virtual int64_t generate();
};

} /* namespace Aux */
#endif /* RANDOMINTEGER_H_ */
