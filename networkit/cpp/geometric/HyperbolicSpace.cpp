/*
 * HyperbolicSpace.cpp
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#include <assert.h>
#include <cmath>


#include "HyperbolicSpace.h"
#include "../auxiliary/Log.h"

using std::abs;
using std::max;

namespace NetworKit {

double HyperbolicSpace::nativeDistance(double firstangle, double firstR, double secondangle, double secondR) {
	assert(firstR >= 0);
	assert(secondR >= 0);
	assert(firstangle >= 0);
	assert(firstangle < 2*M_PI);
	assert(secondangle >= 0);
	assert(secondangle < 2*M_PI);
	double result;
	if (firstangle == secondangle) {
		result = abs(firstR - secondR);
	}
	else {
		double deltaPhi = M_PI - abs(M_PI-abs(firstangle - secondangle));
		double coshDist = cosh(firstR)*cosh(secondR)-sinh(firstR)*sinh(secondR)*cos(deltaPhi);
		if (coshDist >= 1) result = acosh(coshDist);
		else result = 0;
	}
	assert(result >= 0);
	return result;
}

/**
 * This distance measure is taken from the Poincaré disc model.
 */
double HyperbolicSpace::poincareMetric(double phi_a, double  r_a, double phi_b, double r_b) {
	assert(r_a < 1);
	assert(r_b < 1);
	return poincareMetric(polarToCartesian(phi_a, r_a), polarToCartesian(phi_b, r_b));
}

double HyperbolicSpace::poincareMetric(Point2D<double> a, Point2D<double> b) {
	assert(a.length() < 1);
	assert(b.length() < 1);
	double result = acosh( 1 + 2*a.squaredDistance(b) / ((1 - a.squaredLength())*(1 - b.squaredLength())));
	assert(result >= 0);
	return result;
}

//double HyperbolicSpace::nativeHyperbolicDistance(double phi_a, double r_a, double phi_b, double r_b) {
//	/* Returns the hyperbolic distance between points u and v
//	* 2010 paper, eqn: 5
//	*/
//	double deltaPhi = M_PI - abs(M_PI-abs(phi_a - phi_b));
//	double distance = acosh(cosh(r_a)*cosh(r_b) - sinh(r_a)*sinh(r_b)*cos(deltaPhi));
//	return distance;
//}

void HyperbolicSpace::fillPoints(vector<double> &angles, vector<double> &radii, double R, double alpha) {
	fillPoints(angles, radii, 0, 2*M_PI, 0, R, alpha);
}

void HyperbolicSpace::fillPoints(vector<double> &angles, vector<double> &radii, double minPhi, double maxPhi, double minR, double maxR, double alpha) {
	uint64_t n = radii.size();
	assert(angles.size() == n);

	double mincdf = cosh(alpha*minR);
	double maxcdf = cosh(alpha*maxR);
	std::uniform_real_distribution<double> phidist{minPhi, maxPhi};
	std::uniform_real_distribution<double> rdist{mincdf, maxcdf};

	for (uint64_t i = 0; i < n; i++) {
		angles[i] = phidist(Aux::Random::getURNG());
		/**
		 * for the radial coordinate distribution, I took the probability density from Greedy Forwarding in Dynamic Scale-Free Networks Embedded in Hyperbolic Metric Spaces
		 * f (r) = sinh r/(cosh R − 1)
		 * \int sinh = cosh+const
		 */
		double random = rdist(Aux::Random::getURNG());
		radii[i] = (acosh(random)/alpha);
		assert(radii[i] <= maxR);
		assert(radii[i] >= minR);
		assert(angles[i] <= maxPhi);
		assert(angles[i] >= minPhi);
		if (radii[i] == maxR) radii[i] = std::nextafter(radii[i], 0);
		assert(radii[i] < maxR);
	}
}

Point2D<double> HyperbolicSpace::polarToCartesian(double phi, double r) {
	return Point2D<double>(r*cos(phi), r*sin(phi));
}

std::map<index, Point<float> > HyperbolicSpace::polarToCartesian(const vector<double> &angles, const vector<double> &radii) {
	assert(radii.size() == angles.size());
	std::map<index, Point<float> > result;
	for (index i = 0; i < angles.size(); i++) {
		Point2D<double> coord = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
		Point<float> temp(coord[0], coord[1]);
		result.insert(std::make_pair(i, temp));
	}
	return result;
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

void HyperbolicSpace::getEuclideanCircle(Point2D<double> hyperbolicCenter, double hyperbolicRadius, Point2D<double> &euclideanCenter, double &euclideanRadius) {
	double phi_h, r_h;
	HyperbolicSpace::cartesianToPolar(hyperbolicCenter, phi_h, r_h);
	double r_c;
	HyperbolicSpace::getEuclideanCircle(r_h, hyperbolicRadius, r_c, euclideanRadius);
	euclideanCenter = HyperbolicSpace::polarToCartesian(phi_h, r_c);
}

void HyperbolicSpace::getEuclideanCircle(double r_h, double hyperbolicRadius, double &radialCoordOfEuclideanCenter, double &euclideanRadius) {
	double a = cosh(hyperbolicRadius)-1;
	double b = 1-(r_h*r_h);
	radialCoordOfEuclideanCenter = (2*r_h)/(b*a+2);
	euclideanRadius = sqrt(radialCoordOfEuclideanCenter*radialCoordOfEuclideanCenter - (2*r_h*r_h - b*a)/(b*a+2));
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

double HyperbolicSpace::maxRinSlice(double minPhi, double maxPhi, double phi_c, double r_c, double euRadius) {
	double maxCos = max(cos(abs(minPhi - phi_c)), cos(abs(maxPhi - phi_c)));
	//double mirrorAngle;
	//if (phi_c >= M_PI) mirrorAngle = phi_c - M_PI;
	//else mirrorAngle = phi_c + M_PI;

	if (minPhi < phi_c && phi_c < maxPhi) maxCos = 1;
	//applying law of cosines here
	double p = r_c*maxCos;
	double maxR = p + sqrt(p*p - r_c*r_c + euRadius*euRadius);
	return maxR;
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
