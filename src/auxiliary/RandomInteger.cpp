/*
 * RandomInteger.cpp
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "RandomInteger.h"

#include <ctime>

namespace Aux {

std::mt19937_64 RandomInteger::randomEngine = std::mt19937_64(time(nullptr));	// initialize mersenne twister with current time


int64_t RandomInteger::generate(int64_t lower, int64_t upper) {
	std::uniform_int_distribution<int64_t> distribution(lower, upper);
	int64_t i = distribution(RandomInteger::randomEngine);
	return i;
}

int64_t RandomInteger::generateFast(int64_t lower, int64_t upper) {
	int64_t diff = upper - lower + 1;
	int64_t r = rand() % diff;
	r += lower;
	return r;
}

uint64_t RandomInteger::x = 123456789;
uint64_t RandomInteger::y = 362436069;
uint64_t RandomInteger::z = 521288629;

uint64_t RandomInteger::marsaglia() {
	uint64_t t;
	x ^= x << 16;
	x ^= x >> 5;
	x ^= x << 1;
	t = x;
	x = y;
	y = z;
	z = t ^ x ^ y;
	return z;
}

uint64_t RandomInteger::generateFaster(uint64_t upper) {
	return marsaglia() % (upper + 1);
}

uint64_t RandomInteger::generateFaster(uint64_t lower, uint64_t upper) {
	return generateFaster(upper - lower) + lower;
}

} /* namespace Aux */
