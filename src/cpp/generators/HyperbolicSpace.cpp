/*
 * HyperbolicSpace.cpp
 *
 *  Created on: 20.05.2014
 *      Author: moritz
 */

#include <cmath>

#include "HyperbolicSpace.h"

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
	double deltaAngle = abs(firstangle - secondangle); //I don't have to check the direction because of the symmetry of cos
	if (firstR != lastR) {
		lastR = firstR;
		coshlastR = cosh(firstR);
		sinhlastR = sinh(firstR);
	}
	double result = acosh(coshlastR*cosh(secondR) - sinhlastR*sinh(secondR)*cos(deltaAngle));
	return result;
}

double HyperbolicSpace::getDistancePrecached(double firstangle, double firstRcosh, double firstRsinh, double secondangle, double secondRcosh, double secondRsinh) {
	double deltaAngle = abs(firstangle - secondangle); //I don't have to check the direction because of the symmetry of cos
	double result = acosh(firstRcosh*secondRcosh - firstRsinh*secondRsinh*cos(deltaAngle));
	return result;
}
}
