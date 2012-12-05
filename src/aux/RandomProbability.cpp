/*
 * RandomProbability.cpp
 *
 *  Created on: 05.12.2012
 *      Author: cls
 */

#include "RandomProbability.h"

RandomProbability::RandomProbability() {
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	this->uniform = dist;
}

RandomProbability::~RandomProbability() {
	// TODO Auto-generated destructor stub
}

double RandomProbability::generate() {
	double r = this->uniform(this->randomEngine);
	return r;
}
