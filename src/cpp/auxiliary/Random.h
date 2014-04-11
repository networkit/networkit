#ifndef RANDOM_H_
#define RANDOM_H_

/*
 * Random.h
 *
 *  Created on: 02.01.2014
 *      Author: FJW
 */

#include <cstdint>
#include <cassert>
#include <random>

namespace Aux {

/**
 * Provides several functions for random-numbers.
 *
 * All functions are guaranteed to be thread-safe if and only if at least GCC 4.8 is used
 */
namespace Random {

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
 * @returns an integer distributed binomially
 * @param n 	number of trials
 * @param p 	success probability
 */
//uint64_t binomial(double n, double p);

std::mt19937_64& getURNG();

/** @returns a uniform random choice from an indexable container of elements 
*/
template<typename Container> typename Container::value_type& choice(const Container& container) {
    std::uniform_int_distribution<std::size_t> dist{0, container.size() - 1};
    return container[dist(getURNG())];
}


/** @returns a weighted random choice from a vector of elements with given weights */
template <typename Element> Element weightedChoice(const std::vector<std::pair<Element, double> >& weightedElements) {
	double total = 0.0;
	for (uint64_t i = 0; i < weightedElements.size(); i++) {
   		total += weightedElements[i].second;
	}
	double r = Aux::Random::real(total);
	for (uint64_t i = 0; i < weightedElements.size(); i++) {
		if (r < weightedElements[i].second)
	    	return weightedElements[i];
	  		r -= weightedElements[i].second;
	}
	assert(false); // should never get here
}

} // namespace Random
} // namespace Aux


#endif /* RANDOM_H_ */

