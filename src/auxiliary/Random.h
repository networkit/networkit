/*
 * Random.h
 *
 *  Created on: 11.04.2013
 *      Author: cls
 */

#ifndef RANDOM_H_
#define RANDOM_H_

#include <random>

namespace Aux {

class Random {

public:

	Random();

	virtual ~Random();

	/**
	 * @return a uniformly random probability
	 */
	double probability();

	/**
	 * @return a uniformly random integer from the interval [l, u]
	 */
	int64_t integer(int64_t l, int64_t u);

	/**
	 * @return a uniformly random element from a vector
	 */
	template<typename T> T choice(std::vector<T>& vec);

protected:

	std::random_device randomDevice;
	std::default_random_engine randomEngine;

	std::uniform_int_distribution<int64_t> integerDistribution;
	std::uniform_real_distribution<double> probabilityDistribution;


};

} /* namespace Aux */

template<typename T>
inline T Aux::Random::choice(std::vector<T>& vec) {
	int64_t i = this->integer(0, vec.size() - 1);
	return vec[i];
}

#endif /* RANDOM_H_ */
