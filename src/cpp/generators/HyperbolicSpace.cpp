/*
 * HyperbolicSpace.cpp
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#include <assert.h>


#include "HyperbolicSpace.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/Log.h"

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

double HyperbolicSpace::cross(Point<double> a, Point<double> b) {
	return a[0]*b[1] - a[1]*b[0];
}

Point<double> HyperbolicSpace::intersect(Point<double> q, Point<double> s, Point<double> p, Point<double> r) {
	double nominator = cross(q-p, r);
	double denominator = cross (r, s);
	if (denominator == 0) {
		DEBUG("(", r[0], ",", r[1], ") and (", s[0], ",", s[1], ") are colinear or parallel.");
	}
	assert(denominator != 0); //otherwise would be parallel or colinear. Expect better.
	double u = nominator / denominator;
	return q + s.scale(u);
}

Point<double> HyperbolicSpace::circleCenter(Point<double> a, Point<double> b, Point<double> c) {
	Point<double> mid1 = a+b;
	Point<double> mid2 = b+c;
	mid1.scale(0.5);
	mid2.scale(0.5);
	assert(mid1[0] * 2 == a[0] + b[0]);
	assert(mid1[1] * 2 == a[1] + b[1]);
	/**
	 * well, actually this shouldn't happen. If the three points I'm given are in fact only two points,
	 * infinitely many possible circle centers exist.
	 */
	if (a.distance(b) == 0) return mid2;
	if (b.distance(c) == 0) return mid1;
	if (a.distance(c) == 0) return mid1;
	return intersect(mid1, orth(a-b), mid2, orth(b-c));
}

Point<double> HyperbolicSpace::orth(Point<double> a) {
	Point<double> result(-a[1], a[0]);
	return result;
}

Point<double> HyperbolicSpace::mirrorOnCircle(Point<double> a, Point<double> circleCenter,	double radius) {
	double r0 = a.distance(circleCenter);
	double r1 = radius * radius / r0;

	return circleCenter + (a - circleCenter).scale(r1/r0);
}

double HyperbolicSpace::getHyperbolicDistance(Point<double> a, Point<double> b) {
	double phi_1, r_1, phi_2, r_2;
	cartesianToPolar(a, &phi_1, &r_1);
	cartesianToPolar(b, &phi_2, &r_2);
	return getDistance(phi_1, r_1, phi_2, r_2);
}

bool HyperbolicSpace::isBelowArc(Point<double> query, Point<double> a, Point<double> b, double radius) {
	Point<double> origin(0,0);
	Point<double> center = circleCenter(a, b, mirrorOnCircle(a, origin, radius));
	assert(abs((center-a).length() - (center - b).length()) < 0.01);
	return (center-query).length() > (center - a).length();
}

Point<double> HyperbolicSpace::polarToCartesian(double phi, double r) {
	return Point<double>(r*cos(phi), r*sin(phi));
}

void HyperbolicSpace::cartesianToPolar(Point<double> a, double * phi, double *r) {
	*r = a.length();
	if (*r == 0) *phi = 0;
	else if (a[1] >= 0){
		*phi = acos(a[0]/ *r);
	} else {
		*phi = -acos(a[0] / *r);
	}
}
}
