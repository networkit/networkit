/*
 * HyperbolicSpace.cpp
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#include <assert.h>
#include <cmath>


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

double HyperbolicSpace::getHyperbolicDistance(Point2D<double> a, Point2D<double> b) {
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
		(*radii)[i] = hyperbolicRadiusToEuclidean(radius);
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

double HyperbolicSpace::cross(Point2D<double> a, Point2D<double> b) {
	return a[0]*b[1] - a[1]*b[0];
}

Point2D<double> HyperbolicSpace::intersect(Point2D<double> q, Point2D<double> s, Point2D<double> p, Point2D<double> r) {
	double nominator = cross(q-p, r);
	double denominator = cross (r, s);
	if (denominator == 0) {
		DEBUG("(", r[0], ",", r[1], ") and (", s[0], ",", s[1], ") are colinear or parallel.");
	}
	assert(denominator != 0); //otherwise would be parallel or colinear. Expect better.
	double u = nominator / denominator;
	return q + s.scale(u);
}

Point2D<double> HyperbolicSpace::circleCenter(Point2D<double> a, Point2D<double> b, Point2D<double> c) {
	Point2D<double> mid1 = a+b;
	Point2D<double> mid2 = b+c;
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

Point2D<double> HyperbolicSpace::orth(Point2D<double> a) {
	Point2D<double> result(-a[1], a[0]);
	return result;
}

Point2D<double> HyperbolicSpace::mirrorOnCircle(Point2D<double> a, Point2D<double> circleCenter,	double radius) {
	double r0 = a.distance(circleCenter);
	double r1 = radius * radius / r0;

	return circleCenter + (a - circleCenter).scale(r1/r0);
}

bool HyperbolicSpace::isBelowArc(Point2D<double> query, Point2D<double> a, Point2D<double> b, double radius) {
	Point2D<double> origin(0,0);
	Point2D<double> center = circleCenter(a, b, mirrorOnCircle(a, origin, radius));
	assert(abs((center-a).length() - (center - b).length()) < 0.01);
	return (center-query).length() > (center - a).length();
}

Point2D<double> HyperbolicSpace::polarToCartesian(double phi, double r) {
	return Point2D<double>(r*cos(phi), r*sin(phi));
}

void HyperbolicSpace::cartesianToPolar(Point2D<double> a, double &phi, double &r) {
	r = a.length();
	if (r == 0) phi = 0;
	else if (a[1] >= 0){
		phi = acos(a[0]/ r);
	} else {
		phi = -acos(a[0] / r);
	}
	if (phi < 0) phi += 2*M_PI;
}

void HyperbolicSpace::getTransmutationCircle(Point2D<double> source,
		Point2D<double> target, double R, Point2D<double> &circleCenter, double &circleRadius) {
	/**
	 * make sure target is outwards
	 */
	if (source.length() > target.length()) {
		Point2D<double> temp = target;
		target = source;
		source = temp;
	}

	double dist = target.distance(source);
	double lambdanom = (-(source[0]*source[0]) - (source[1]*source[1]) + R*R);
	double lambdadenom = dist * dist+ 2*(source[0]*(target[0] - source[0])+source[1]*(target[1]-source[1]));

	circleCenter = (target - source).scale(lambdanom/lambdadenom) + source;

	circleRadius = pow(target.distance(circleCenter) * source.distance(circleCenter), 0.5);
}

double HyperbolicSpace::hyperbolicDistanceToArc(Point2D<double> query,
		Point2D<double> a, Point2D<double> b, double R) {

	assert(R == 1);

	/**
	 * get direct distances
	 */

	double qToA = getHyperbolicDistance(query,a);
	double qToB = getHyperbolicDistance(query,b);

	/**
	 * transform
	 */

	Point2D<double> origin(0,0);
	Point2D<double> m;
	double radius;
	getTransmutationCircle(origin, query, R, m, radius);

	Point2D<double> adash = mirrorOnCircle(a, m, radius);
	Point2D<double> bdash = mirrorOnCircle(b,m,radius);
	Point2D<double> querydash = mirrorOnCircle(query,m,radius);
	assert(querydash.length() < 0.0001);//should be origin

	/**
	* get connecting arc
	 */

	Point2D<double> c = circleCenter(adash, bdash, mirrorOnCircle(adash, origin, R));
	double arcRadius = c.distance(adash);
	double phi_c, r_c;
	cartesianToPolar(c, phi_c, r_c);
	double directDistance = r_c - arcRadius;

	Point2D<double> closestOnCircle = polarToCartesian(phi_c, r_c - arcRadius);
	assert(closestOnCircle.distance(c) - arcRadius < 0.0001);
	double aToCoC = getHyperbolicDistance(adash, closestOnCircle);
	double bToCoC = getHyperbolicDistance(bdash, closestOnCircle);
	double aToB = getHyperbolicDistance(adash, bdash);

	if (aToCoC < aToB && bToCoC < aToB) return directDistance;
	else return std::min(qToA, qToB);
}

Point2D<double> HyperbolicSpace::getPointOnHyperbolicCircle(Point2D<double> hyperbolicCenter, double radius) {
	double phi_q, r_q;
	HyperbolicSpace::cartesianToPolar(hyperbolicCenter, phi_q, r_q);
	double hyperbolicRadius = HyperbolicSpace::EuclideanRadiusToHyperbolic(r_q);
	if (hyperbolicRadius < radius) {
		//origin is part of circle, need to use other method
		//use right triangle at origin to construct new point at phi_p, r_p
		double phi_p, r_p;
		phi_p = phi_q + M_PI/2;
		if (phi_p > 2*M_PI) phi_p -= 2*M_PI;
		double r_p_nom = (cosh(radius)-1)*(1-r_q*r_q)-2*r_q*r_q;
		double r_p_denom = (2+(cosh(radius)-1)*(1-r_q*r_q));
		r_p = pow(r_p_nom/r_p_denom, 0.5);
		return HyperbolicSpace::polarToCartesian(phi_p, r_p);
	} else {
		double sqL = hyperbolicCenter.squaredLength();
		double nom = (cosh(radius)-1)*(1-sqL)*(1-sqL);
		double denom = 4*sqL;
		double gamma = acos(1-nom/denom);
		double newphi;
		if (gamma + phi_q < 2*M_PI) {
			newphi = gamma + phi_q;
		} else {
			assert(phi_q - gamma >= 0);
			newphi = phi_q - gamma;
		}
		return HyperbolicSpace::polarToCartesian(newphi, r_q);
	}
}

void HyperbolicSpace::getEuclideanCircle(Point2D<double> hyperbolicCenter, Point2D<double> pointOnEdge, Point2D<double> &euclideanCenter, double &euclideanRadius) {
	Point2D<double> origin;
	Point2D<double> center = circleCenter(hyperbolicCenter, pointOnEdge, mirrorOnCircle(hyperbolicCenter, origin, 1));
	Point2D<double> edgeCenterVector = center - pointOnEdge;
	Point2D<double> arc = hyperbolicCenter - origin;
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

double HyperbolicSpace::hyperbolicSpaceInEuclideanCircle(double r_c, double d_c,
		double r_max) {
	double result = 0;
	assert(r_c >= 0);
	assert(d_c >= 0);
	assert(r_c <= r_max);
	double min = r_c - d_c;
	double max = std::min(r_c+d_c, r_max);

	if (d_c > r_c) {
		//the query circle overlaps the origin

		if (d_c - r_c < r_max) {
			//remaining query circle is smaller than the disk representation
			result += 2*M_PI*(cosh(EuclideanRadiusToHyperbolic(d_c-r_c))-1);//adding small circle around origin
		} else {
			result += 2*M_PI*(cosh(EuclideanRadiusToHyperbolic(r_max))-1);//adding small circle around origin
		}
		assert(result <= 2*M_PI*(cosh(EuclideanRadiusToHyperbolic(r_max))-1));
		min = std::nextafter(d_c-r_c, std::numeric_limits<double>::max());//correcting integral start to exclude circle
	}

	/**
	 * Now, the integral.
	 * It is 4\int_{min}^{max} \text{acos}(\frac{r_c^2-d_c^2+r^2}{2r_c\cdot r})  \cdot \frac{1}{1-r^2} \cdot (\sinh (\text{acosh}( 1 + 2\frac{r^2}{1 - r^2})))\,dr
	 * This solution was computed by WolframAlpha
	 */

	if (max < min) return result;

	auto realpart = [](double r, double d, double c) {
		double result = acos((c*c-d*d+r*r) / (2*c*r)) / (r*r-1);
		return result;
	};

	/**
	 * Maybe the denominator can be omitted, atan2 is probably the same.
	 */
	auto firstlogpart = [](double r, double d, double c) {
		double s = (c*c-d*d);
		//double denominator = r*r*s*s;
		double rsqs = r*r+s;
		double real = -2*s*sqrt(4*c*c*r*r-rsqs*rsqs);
		double imag = -4*c*c*r*r+2*s*r*r+2*s*s;
		return atan2(imag, real)/2;
	};

	auto secondlogpart = [](double r, double d, double c) {
		double s = (c*c-d*d);
		double rsqs = r*r+s;
		//double denominator = (r*r-1)*(s-1);
		double real = sqrt(4*c*c*r*r-rsqs*rsqs);
		double imag = 2*c*c*(r*r+1)-(s+1)*rsqs;
		imag = imag / sqrt((s+1)*(s+1)-(4*c*c));
		return (s-1)*atan2(imag, real)/(2*sqrt((s+1)*(s+1)-(4*c*c)));
	};

	double lower = -realpart(min, d_c, r_c) -firstlogpart(min, d_c, r_c) + secondlogpart(min, d_c, r_c);
	double upper = -realpart(max, d_c, r_c) -firstlogpart(max, d_c, r_c) + secondlogpart(max, d_c, r_c);
	result += 4*(upper - lower);
	return result;
}
}
