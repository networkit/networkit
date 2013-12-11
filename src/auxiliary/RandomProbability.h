/*
 * RandomProbability.h
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef RANDOMPROBABILITY_H_
#define RANDOMPROBABILITY_H_

#include <random>

namespace Aux {

/**
 * Random probability value generators.
 */
class RandomProbability {

private:


	static std::mt19937_64 randomEngine; 
	static uint64_t x, y, z;	// Marsaglia's generator needs these states

	RandomProbability() = delete;


public:

	static double generate();

	/**
	 * @return Random probability in half-open interval [0 ... 1)
	 */
	static inline double generateFast() {
		return ((double) rand()) / (RAND_MAX + 1.0);
	}

	/**
	 * @return Random float in [0 ... 1)
	 */
	static float randomFloat();
};

}

#endif /* RANDOMPROBABILITY_H_ */
