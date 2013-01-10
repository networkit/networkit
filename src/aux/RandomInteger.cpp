/*
 * RandomInteger.cpp
 *
 *  Created on: 10.01.2013
 *      Author: cls
 */

#include "RandomInteger.h"

namespace EnsembleClustering {

RandomInteger::RandomInteger(int64_t lower, int64_t upper) : randomEngine((unsigned int) time(0)), distribution(lower, upper) {
}

RandomInteger::~RandomInteger() {
	// TODO Auto-generated destructor stub
}

int64_t RandomInteger::generate() {
	int64_t i = this->distribution(this->randomEngine);
	return i;
}

} /* namespace EnsembleClustering */
