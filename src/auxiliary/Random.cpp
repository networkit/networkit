/*
 * Random.cpp
 *
 *  Created on: 02.01.2014
 *      Author: FJW
 */

#include "Random.h"

#include <random>
#include <cmath>

namespace Aux {
namespace Random {

namespace {
std::mt19937_64& getURNG();
}

uint64_t getSeed() {
	thread_local static std::random_device urng{};
	std::uniform_int_distribution<uint64_t> dist{};
	return dist(urng);
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

namespace {
std::mt19937_64& getURNG() {
	thread_local static std::mt19937_64 generator{getSeed()};
	return generator;
}
} // anonymous namespace

} // namespace Random
} // namespace Aux
