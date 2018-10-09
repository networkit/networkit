/*
 * Random.cpp
 *
 *  Created on: 02.01.2014
 *      Author: FJW
 */

#include <cmath>
#include <omp.h>
#include <limits>
#include <atomic>

#include "Random.h"

namespace Aux {
namespace Random {

static std::atomic<uint64_t> seedValue{0};
static std::atomic<bool> seedUseThreadId{false};

// The global seed generation functions as an epoch counter
// for each time, the seed was updated. At the same time it
// releases the changes made to the seed values to other
// threads which always access it by an aquire.
// seedValues and seedUseThreadId can hence be accesses
// relaxed.
// As long as globalSeedGeneration is zero, the getURNG()
// ignores these two variables.
static std::atomic<uint64_t> globalSeedGeneration{0}; //< global seed generation, updated on every setSeed-call

void setSeed(uint64_t seed, bool useThreadId) {
	seedValue.store(seed, std::memory_order_relaxed);
	seedUseThreadId.store(useThreadId, std::memory_order_relaxed);
	globalSeedGeneration.fetch_add(1, std::memory_order_release) ;
	getURNG(); // update local seed value
}

static uint64_t getSeed_(uint64_t globSeedGen) {
	if (!globSeedGen) {
		thread_local static std::random_device urng{};
		std::uniform_int_distribution<uint64_t> dist{};
		return dist(urng);

	} else if (seedUseThreadId.load(std::memory_order_relaxed)) {
		return seedValue.load(std::memory_order_relaxed) + omp_get_thread_num();

	} else {
		return seedValue.load(std::memory_order_relaxed);
	}
}

uint64_t getSeed() {
	return getSeed_(globalSeedGeneration.load(std::memory_order_acquire));
}

bool getUseThreadId() {
	return seedUseThreadId;
}

std::mt19937_64& getURNG() {
	thread_local static std::mt19937_64 generator{getSeed()};
	thread_local static uint64_t localSeedGeneration{0};

	auto globSeedGen = globalSeedGeneration.load(std::memory_order_acquire);

	if (localSeedGeneration != globSeedGen) {
		generator.seed(getSeed_(globSeedGen));
		localSeedGeneration = globSeedGen;
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

} // namespace Random
} // namespace Aux
