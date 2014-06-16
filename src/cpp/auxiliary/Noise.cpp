/*
 * Noise.cpp
 *
 *  Created on: 29.11.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "Noise.h"

#include "Random.h"

namespace Aux {


Noise::Noise(double l, double u) : uniform{l, u} {}

double Noise::add(double x) {
	double r = this->uniform(Random::getURNG());
	return x + r;
}

}
