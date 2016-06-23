#ifndef RANDOM_H_
#define RANDOM_H_

/*
 * Random.h
 *
 *  Created on: 02.01.2014
 *      Author: FJW
 */

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include "Log.h"

namespace Aux {

/**
 * Provides several functions for random-numbers.
 *
 * All functions are guaranteed to be thread-safe if and only if at least GCC 4.8 is used
 */
namespace Random {

/**
 * Sets the random seed that is used globally.
 * @param seed The seed value
 * @param useThreadId If the thread id shall be added to the seed value
 */
void setSeed(uint64_t seed, bool useThreadId);

/**
 * @returns a high-quality random seed for an URNG.
 */
uint64_t getSeed();

/**
 * @returns a reference to a seeded URNG that is thread_local iff GCC 4.8 or later is used.
 */
std::mt19937_64& getURNG();


/**
 * @returns an integer distributed uniformly in an inclusive range;
 * @param upperBound the upper bound, default = UNINT64_T_MAX
 * @param lowerBound the lower bound, default = 0
 */
uint64_t integer();
uint64_t integer(uint64_t upperBound);
uint64_t integer(uint64_t lowerBound, uint64_t upperBound);

/**
 * @returns a double distributed uniformly in a half-open range: [lowerBound, upperBound)
 * @param upperBound default = 1.0
 * @param lowerBound default = 0.0
 */
double real();
double real(double upperBound);
double real(double lowerBound, double upperBound);

/**
 * @returns a double distributed uniformly in the range [0, 1]
 * @note this differs from real() in returning a value in a closed instead of a half-open range.
 */
double probability();

/**
 * @returns a size_t in the range [0, max - 1] that can for example be used to
 * access a random element of a container.
 */
std::size_t index(std::size_t max);

/**
 * @returns an integer distributed binomially
 * @param n 	number of trials
 * @param p 	success probability
 */
//uint64_t binomial(double n, double p);


/**
 * @returns a uniform random choice from an indexable container of elements.
 */
template<typename Container>
typename Container::const_reference choice(const Container& container) {
	return container[index(container.size())];
}


/**
 * @returns a weighted random choice from a vector of elements with given weights.
 */
template <typename Element>
const Element& weightedChoice(const std::vector<std::pair<Element, double>>& weightedElements) {
	if (weightedElements.size() == 0)
		throw std::runtime_error("Random::weightedChoice: input size equal to 0");
	double total = 0.0;
	for (const auto& entry : weightedElements) {
		assert(entry.second >= 0.0 && "This algorithm only works with non-negative weights");
		total += entry.second;
	}
	double r = Aux::Random::real(total);
	for (const auto& entry : weightedElements) {
		if (r < entry.second) {
			return entry.first;
		}
		r -= entry.second;
	}
	throw std::runtime_error("Random::weightedChoice: should never get here"); // should never get here
}

} // namespace Random
} // namespace Aux


#endif /* RANDOM_H_ */
