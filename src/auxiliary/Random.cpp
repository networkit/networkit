/*
 * Random.cpp
 *
 *  Created on: 11.04.2013
 *      Author: cls
 */

#include "Random.h"

namespace Aux {

Random::Random() : randomEngine(this->randomDevice()), probabilityDistribution(0.0, 1.0), integerDistribution() {
	// TODO Auto-generated constructor stub

}

Random::~Random() {
	// TODO Auto-generated destructor stub
}

double Random::probability() {
	double r = this->probabilityDistribution(this->randomEngine);
	return r;
}

int64_t Random::integer(int64_t l, int64_t u) {
	int64_t i = this->integerDistribution(this->randomEngine, std::uniform_int_distribution<int64_t>{l, u}.param());
	return i;
}

} /* namespace Aux */
