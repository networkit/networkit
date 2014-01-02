/*
 * Noise.cpp
 *
 *  Created on: 29.11.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "Noise.h"

namespace Aux {


Noise::Noise(double l, double u) :
	uniform{l, u}, randomEngine{std::random_device{}()} {}

double Noise::add(double x) {
	double r = this->uniform(this->randomEngine);
	return x + r;
}

}
