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

class RandomProbability {

protected:

	std::random_device randomDevice;
	std::default_random_engine randomEngine;
	std::uniform_real_distribution<double> distribution;


public:

	RandomProbability();

	virtual ~RandomProbability();

	double generate();

	/**
	 * @return Random probability in half-open interval [0 ... 1)
	 */
	inline double generateFast() {
		return ((double) rand()) / (RAND_MAX + 1.0);
	}

	/**
	 * @return Random float in [0 ... 1)
	 */
	float randomFloat() const;
};

}

#endif /* RANDOMPROBABILITY_H_ */
