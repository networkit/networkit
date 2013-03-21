/*
 * RandomInteger.cpp
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "RandomInteger.h"

namespace Aux {

RandomInteger::RandomInteger(int64_t lower, int64_t upper) : randomEngine(this->randomDevice()), distribution(lower, upper) {
}

RandomInteger::~RandomInteger() {
	// TODO Auto-generated destructor stub
}

int64_t RandomInteger::generate() {
	int64_t i = this->distribution(this->randomEngine);
	return i;
}

} /* namespace Aux */
