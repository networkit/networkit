/*
 * HyperbolicSpace.cpp
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#include <cmath>
#include <assert.h>

#include "HyperbolicSpace.h"
#include "../auxiliary/Random.h"

using std::abs;

namespace NetworKit {

double HyperbolicSpace::lastR = 0;
double HyperbolicSpace::coshlastR = 1;
double HyperbolicSpace::sinhlastR = 0;

HyperbolicSpace::HyperbolicSpace() {
	// TODO Auto-generated constructor stub
	radius = 1;
}

HyperbolicSpace::~HyperbolicSpace() {
	// TODO Auto-generated destructor stub
}

HyperbolicSpace::HyperbolicSpace(double R) {
	radius = R;
}

double HyperbolicSpace::getRadius() {
	return radius;
}

/**
 * This distance measure is taken from the Poincaré disc model.
 */
double HyperbolicSpace::getDistance(double firstangle, double firstR, double secondangle, double secondR) {
	if (firstangle == secondangle && firstR == secondR) return 0;
	double deltaAngle = abs(firstangle - secondangle); //I don't have to check the direction because of the symmetry of cos
	if (firstR != lastR) {
		lastR = firstR;
		coshlastR = cosh(firstR);
		sinhlastR = sinh(firstR);
	}
	double result = acosh(coshlastR*cosh(secondR) - sinhlastR*sinh(secondR)*cos(deltaAngle));
	assert(result >= 0);
	return result;
}

//TODO: add const keywords where appropriate
void HyperbolicSpace::fillPoints(vector<double> * angles, vector<double> * radii, double stretch, double alpha) {
	uint64_t n = radii->size();
	double R = stretch*acosh((double)n/(2*M_PI)+1);
	assert(angles->size() == n);
	for (uint64_t i = 0; i < n; i++) {
		(*angles)[i] = Aux::Random::real(0, 2*M_PI);
		/**
		 * for the radial coordinate distribution, I took the probability density from Greedy Forwarding in Dynamic Scale-Free Networks Embedded in Hyperbolic Metric Spaces
		 * f (r) = sinh r/(cosh R − 1)
		 * \int sinh = cosh+const
		 */
		//TODO: implement alpha
		double maxcdf = cosh(alpha*R);
		double random = Aux::Random::real(1, maxcdf);
		(*radii)[i] = acosh(random)/alpha;
		assert((*radii)[i] <= R);
	}
}

double HyperbolicSpace::getDistancePrecached(double firstangle, double firstRcosh, double firstRsinh, double secondangle, double secondRcosh, double secondRsinh) {
	double deltaAngle = abs(firstangle - secondangle); //I don't have to check the direction because of the symmetry of cos
	if (deltaAngle == 0 && firstRcosh == secondRcosh && firstRsinh == secondRsinh) return 0;//points are identical, sometimes resulted in -nan
	double result = acosh(firstRcosh*secondRcosh - firstRsinh*secondRsinh*cos(deltaAngle));
	assert(result >= 0);
	return result;
}
}
