/*
 * Random.cpp
 *
 *  Created on: 02.01.2014
 *      Author: FJW
 */

#include <cmath>
#include <random>

#include "Random.h"

// If GCC does not support thread local, we are sad and don't use it:
#ifdef __GNUC__
#	if (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR >= 8)
#		define AUX_THREAD_LOCAL thread_local
#	else
#		define AUX_THREAD_LOCAL
#	endif
#else // we don't know our plattform, so don't support it:
#	define AUX_THREAD_LOCAL 
#endif

namespace Aux {
namespace Random {

namespace {

/**
 * private helper that provides access to a thread-local URNG (uniform random number generator)
 */
std::mt19937_64& getURNG();
}

uint64_t getSeed() {
	AUX_THREAD_LOCAL static std::random_device urng{};
	std::uniform_int_distribution<uint64_t> dist{};
	return dist(urng);
}

uint64_t integer() {
	AUX_THREAD_LOCAL static std::uniform_int_distribution<uint64_t> dist{};
	return dist(getURNG());
}
uint64_t integer(uint64_t upperBound) {
	std::uniform_int_distribution<uint64_t> dist{0, upperBound};
	return dist(getURNG());
}
uint64_t integer(uint64_t lowerBound, uint64_t upperBound) {
	std::uniform_int_distribution<uint64_t> dist{lowerBound, upperBound};
	return dist(getURNG());
}

double real() {
	AUX_THREAD_LOCAL static std::uniform_real_distribution<double> dist{};
	return dist(getURNG());
}
double real(double upperBound) {
	std::uniform_real_distribution<double> dist{0.0, upperBound};
	return dist(getURNG());
}
double real(double lowerBound, double upperBound) {
	std::uniform_real_distribution<double> dist{lowerBound, upperBound};
	return dist(getURNG());
}

double probability() {
	AUX_THREAD_LOCAL static std::uniform_real_distribution<double> dist{0.0, std::nexttoward(1.0, 2.0)};
	return dist(getURNG());
}

uint64_t binomial(uint64_t n, double p) {
	AUX_THREAD_LOCAL static std::binomial_distribution<uint64_t> dist{n, p};
	return dist(getURNG());
}

namespace {
std::mt19937_64& getURNG() {
	AUX_THREAD_LOCAL static std::mt19937_64 generator{getSeed()};
	return generator;
}
} // anonymous namespace

} // namespace Random
} // namespace Aux
