/*
 * Noise.cpp
 *
 *  Created on: 29.11.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "Noise.h"

namespace Aux {


Noise::Noise(double l, double u) {
	this->lowerBound = l;
	this->upperBound = u;
	std::uniform_real_distribution<double> dist(l, u);
	this->uniform = dist;
}

Noise::~Noise() {
	// TODO Auto-generated destructor stub
}

double Noise::add(double x) {
	double r = this->uniform(this->randomEngine);
	return x + r;
}

}
