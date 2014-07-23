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
using std::max;

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
double HyperbolicSpace::getHyperbolicDistance(double phi_a, double  r_a, double phi_b, double r_b) {
	/**
	 * quick and dirty to see if it works. TODO: clean up later.
	 */
	assert(r_a < 1);
	assert(r_b < 1);
	return getHyperbolicDistance(polarToCartesian(phi_a, r_a), polarToCartesian(phi_b, r_b));
}

double HyperbolicSpace::getHyperbolicDistance(Point<double> a, Point<double> b) {
	double result = acosh( 1 + 2*a.squaredDistance(b) / ((1 - a.squaredLength())*(1 - b.squaredLength())));
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
		double radius = (acosh(random)/alpha);
		//now translate into coordinates of Poincaré disc
		double rad_nom = (cosh(radius)-1);
		double rad_denom = (cosh(radius)+1);
		(*radii)[i] = sqrt(rad_nom/rad_denom);
		assert((*radii)[i] < 1);
	}
}
/**
double HyperbolicSpace::getDistancePrecached(double firstangle, double firstRcosh, double firstRsinh, double secondangle, double secondRcosh, double secondRsinh) {
	double deltaAngle = abs(firstangle - secondangle); //I don't have to check the direction because of the symmetry of cos
	if (deltaAngle == 0 && firstRcosh == secondRcosh && firstRsinh == secondRsinh) return 0;//points are identical, sometimes resulted in -nan
	double result = acosh(firstRcosh*secondRcosh - firstRsinh*secondRsinh*cos(deltaAngle));
	assert(result >= 0);
	return result;
}
*/

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

bool HyperbolicSpace::isBelowArc(Point<double> query, Point<double> a, Point<double> b, double radius) {
	Point<double> origin(0,0);
	Point<double> center = circleCenter(a, b, mirrorOnCircle(a, origin, radius));
	assert(abs((center-a).length() - (center - b).length()) < 0.01);
	return (center-query).length() > (center - a).length();
}

Point<double> HyperbolicSpace::polarToCartesian(double phi, double r) {
	return Point<double>(r*cos(phi), r*sin(phi));
}

void HyperbolicSpace::cartesianToPolar(Point<double> a, double &phi, double &r) {
	r = a.length();
	if (r == 0) phi = 0;
	else if (a[1] >= 0){
		phi = acos(a[0]/ r);
	} else {
		phi = -acos(a[0] / r);
	}
	if (phi < 0) phi += 2*M_PI;
}

void HyperbolicSpace::getTransmutationCircle(Point<double> source,
		Point<double> target, double R, Point<double> &circleCenter, double &circleRadius) {
	/**
	 * make sure target is outwards
	 */
	if (source.length() > target.length()) {
		Point<double> temp = target;
		target = source;
		source = temp;
	}

	double dist = target.distance(source);
	double lambdanom = (-(source[0]*source[0]) - (source[1]*source[1]) + R*R);
	double lambdadenom = dist * dist+ 2*(source[0]*(target[0] - source[0])+source[1]*(target[1]-source[1]));

	circleCenter = (target - source).scale(lambdanom/lambdadenom) + source;

	circleRadius = pow(target.distance(circleCenter) * source.distance(circleCenter), 0.5);
}

double HyperbolicSpace::hyperbolicDistanceToArc(Point<double> query,
		Point<double> a, Point<double> b, double R) {

	assert(R == 1);

	/**
	 * get direct distances
	 */

	double qToA = getHyperbolicDistance(query,a);
	double qToB = getHyperbolicDistance(query,b);

	/**
	 * transform
	 */

	Point<double> origin(0,0);
	Point<double> m;
	double radius;
	getTransmutationCircle(origin, query, R, m, radius);

	Point<double> adash = mirrorOnCircle(a, m, radius);
	Point<double> bdash = mirrorOnCircle(b,m,radius);
	Point<double> querydash = mirrorOnCircle(query,m,radius);
	assert(querydash.length() < 0.0001);//should be origin

	/**
	* get connecting arc
	 */

	Point<double> c = circleCenter(adash, bdash, mirrorOnCircle(adash, origin, R));
	double arcRadius = c.distance(adash);
	double phi_c, r_c;
	cartesianToPolar(c, phi_c, r_c);
	double directDistance = r_c - arcRadius;

	Point<double> closestOnCircle = polarToCartesian(phi_c, r_c - arcRadius);
	assert(closestOnCircle.distance(c) - arcRadius < 0.0001);
	double aToCoC = getHyperbolicDistance(adash, closestOnCircle);
	double bToCoC = getHyperbolicDistance(bdash, closestOnCircle);
	double aToB = getHyperbolicDistance(adash, bdash);

	if (aToCoC < aToB && bToCoC < aToB) return directDistance;
	else return std::min(qToA, qToB);
}

Point<double> HyperbolicSpace::getPointOnHyperbolicCircle(Point<double> hyperbolicCenter, double radius) {
	double sqL = hyperbolicCenter.squaredLength();
	double nom = (cosh(radius)-1)*(1-sqL)*(1-sqL);
	double denom = 4*sqL;
	double gamma = acos(1-nom/denom);
	double phi, r, newphi;
	HyperbolicSpace::cartesianToPolar(hyperbolicCenter, phi, r);
	if (gamma + phi < 2*M_PI) {
		newphi = gamma + phi;
	} else {
		assert(phi - gamma >= 0);
		newphi = phi - gamma;
	}
	return HyperbolicSpace::polarToCartesian(newphi, r);
}

void HyperbolicSpace::getEuclideanCircle(Point<double> hyperbolicCenter, Point<double> pointOnEdge, Point<double> &euclideanCenter, double &euclideanRadius) {
	Point<double> origin;
	Point<double> center = circleCenter(hyperbolicCenter, pointOnEdge, mirrorOnCircle(hyperbolicCenter, origin, 1));
	Point<double> edgeCenterVector = center - pointOnEdge;
	Point<double> arc = hyperbolicCenter - origin;
	euclideanCenter = intersect(pointOnEdge, orth(edgeCenterVector), origin, arc);
	euclideanRadius = euclideanCenter.distance(pointOnEdge);
}

double HyperbolicSpace::hyperbolicRadiusToEuclidean(double hyperbolicRadius) {
	double ch = cosh(hyperbolicRadius);
	return sqrt((ch-1)/(ch+1));
}

double HyperbolicSpace::EuclideanRadiusToHyperbolic(double euclideanRadius) {
	double eusq = euclideanRadius*euclideanRadius;
	double result = acosh( 1 + 2*eusq / ((1 - eusq)));
	return result;
}
}
