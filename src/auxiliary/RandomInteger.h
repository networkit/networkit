/*
 * RandomInteger.h
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef RANDOMINTEGER_H_
#define RANDOMINTEGER_H_

#include <random>

namespace Aux {

class RandomInteger {



private:

	static std::mt19937_64 randomEngine; 
	static uint64_t x, y, z;	// Marsaglia's generator needs these states

	RandomInteger() = delete;

	/**
	* Marsaglia's xorshf generator, see: http://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c/1640399#1640399
	*/
	static uint64_t marsaglia();

public:

	/**
	 * Generate random integer in closed interval [lower, upper].
	 */
	static int64_t generate(int64_t lower, int64_t upper);

	/**
	 * Generate random integer in closed interval [lower, upper]. Faster generation,
	 * but without large strength in terms of randomness;
	 */
	static int64_t generateFast(int64_t lower, int64_t upper);


	static uint64_t generateFaster(uint64_t upper);

	/**
	 * EXPERIMENTAL: Generate random integer in closed interval [lower, upper].
	 */
	static uint64_t generateFaster(uint64_t lower, uint64_t upper);
};

} /* namespace Aux */
#endif /* RANDOMINTEGER_H_ */
