/*
 * Noise.h
 *
 *  Created on: 29.11.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOISE_H_
#define NOISE_H_

#include <random>

namespace Aux {

/**
 * Noise is random addition to a signal. This class provides methods
 * which add random numbers to their inputs in order to enable randomization.
 */
class Noise {
protected:
	std::uniform_real_distribution<double> uniform;
	
public:
	/**
	 * @param[in]	l	lower bound for added random number
	 * @param[in]	u	upper bound for added random number
	 */
	Noise(double l, double u);
	
	virtual ~Noise() = default;
	
	/**
	 * Add noise to double.
	 *
	 * @param[in]	 x	input
	 * @return	input plus noise
	 */
	double add(double x);
	
	
};


}
#endif /* NOISE_H_ */
