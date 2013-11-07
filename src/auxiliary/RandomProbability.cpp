/*
 * RandomProbability.cpp
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "RandomProbability.h"

namespace Aux {

std::mt19937_64 RandomProbability::randomEngine = std::mt19937_64(time(nullptr));	// initialize mersenne twister with current time



double RandomProbability::generate() {
	std::uniform_real_distribution<double> distribution;
	double r = distribution(RandomProbability::randomEngine);
	return r;
}

float RandomProbability::randomFloat() {
	return ((float) rand()) / (RAND_MAX + 1.0f);
}


}
