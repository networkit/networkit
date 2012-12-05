/*
 * RandomProbability.h
 *
 *  Created on: 05.12.2012
 *      Author: cls
 */

#ifndef RANDOMPROBABILITY_H_
#define RANDOMPROBABILITY_H_

#include <random>

class RandomProbability {

protected:

	std::uniform_real_distribution<double> uniform;
	std::default_random_engine randomEngine;

public:

	RandomProbability();

	virtual ~RandomProbability();

	virtual double generate();
};

#endif /* RANDOMPROBABILITY_H_ */
