/*
 * RandomInteger.cpp
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "RandomInteger.h"

namespace Aux {

RandomInteger::RandomInteger() : randomEngine(this->randomDevice()), distribution() {
}

RandomInteger::~RandomInteger() {
	// TODO Auto-generated destructor stub
}

int64_t RandomInteger::generate(int64_t lower, int64_t upper) {
	int64_t i = this->distribution(this->randomEngine, std::uniform_int_distribution<int64_t>{lower, upper}.param());
	return i;
}

} /* namespace Aux */
