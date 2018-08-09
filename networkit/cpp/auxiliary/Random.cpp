/*
 * Random.cpp
 *
 *  Created on: 02.01.2014
 *      Author: FJW
 */

#include <cmath>
#include <omp.h>
#include <limits>

#include "Random.h"

namespace Aux {
namespace Random {

static bool staticSeed = false;
static uint64_t seedValue = 0;
static uint64_t globalSeedGeneration = 0; // global seed generation, updated on every setSeed-call
static bool seedUseThredId = false;

void setSeed(uint64_t seed, bool useThreadId) {
	seedValue = seed;
	staticSeed = true;
	seedUseThredId = useThreadId;
	++globalSeedGeneration;
	getURNG(); // update local seed value
}

uint64_t getSeed() {
	if (!staticSeed) {
		thread_local static std::random_device urng{};
		std::uniform_int_distribution<uint64_t> dist{};
		return dist(urng);
	} else if (seedUseThredId) {
		return seedValue + omp_get_thread_num();
	} else {
		return seedValue;
	}
}

std::mt19937_64& getURNG() {
	thread_local static std::mt19937_64 generator{getSeed()};
	thread_local static uint64_t localSeedGeneration = std::numeric_limits<uint64_t>::max();
	if (staticSeed && localSeedGeneration != globalSeedGeneration) {
		generator.seed(getSeed());
		localSeedGeneration = globalSeedGeneration;
	}
	return generator;
}

uint64_t integer() {
	thread_local static std::uniform_int_distribution<uint64_t> dist{};
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
	thread_local static std::uniform_real_distribution<double> dist{};
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
	thread_local static std::uniform_real_distribution<double> dist{0.0, std::nexttoward(1.0, 2.0)};
	return dist(getURNG());
}

std::size_t index(std::size_t max) {
	assert(max > 0 && "There have to be valid indexes");
	std::uniform_int_distribution<std::size_t> dist{0, max - 1};
	return dist(getURNG());
}

// uint64_t binomial(double n, double p) {
// 	std::binomial_distribution<uint64_t> dist(n, p);
// 	return dist(getURNG());
// }

} // namespace Random
} // namespace Aux
