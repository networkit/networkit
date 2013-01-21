/*
 * RandomProbability.cpp
 *
 *  Created on: 05.12.2012
 *      Author: cls
 */

#include "RandomProbability.h"

namespace Aux {

RandomProbability::RandomProbability() :  randomEngine((unsigned int) time(0)), distribution(0.0, 1.0) {
}

RandomProbability::~RandomProbability() {
	// TODO Auto-generated destructor stub
}

double RandomProbability::generate() {
	double r = this->distribution(this->randomEngine);
	return r;
}

}
