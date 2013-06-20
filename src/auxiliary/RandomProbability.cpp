/*
 * RandomProbability.cpp
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "RandomProbability.h"

namespace Aux {

RandomProbability::RandomProbability() :  randomEngine(this->randomDevice()), distribution(0.0, 1.0) {
}

RandomProbability::~RandomProbability() {
	// TODO Auto-generated destructor stub
}

double RandomProbability::generate() {
	double r = this->distribution(this->randomEngine);
	return r;
}

double RandomProbability::generateFast() {
	return ((double) rand()) / (RAND_MAX + 1.0);
}

// TODO: for faster generation, use rand() function

float RandomProbability::randomFloat() const {
	return ((float) rand()) / (RAND_MAX + 1.0f);
}


}
