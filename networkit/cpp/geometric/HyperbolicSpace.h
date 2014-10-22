/*
 * HyperbolicSpace.h
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef HYPERBOLICSPACE_H_
#define HYPERBOLICSPACE_H_

#include <vector>
#include <cmath>
#include "Point2D.h"

using std::vector;
using std::abs;

namespace NetworKit {

class HyperbolicSpace {
public:
	HyperbolicSpace();
	virtual ~HyperbolicSpace();
	static void fillPoints(vector<double> * angles, vector<double> * radii, double stretch, double alpha);
	static double getHyperbolicDistance(double firstangle, double firstR, double secondangle, double secondR);
	static double getHyperbolicDistance(Point2D<double> a, Point2D<double> b);
	static Point2D<double> polarToCartesian(double phi, double r);
	static void cartesianToPolar(Point2D<double> a, double &phi, double &r);
	static void getEuclideanCircle(Point2D<double> hyperbolicCenter, double hyperbolicRadius, Point2D<double> &euclideanCenter, double &euclideanRadius);
	static void getEuclideanCircle(double r_h, double hyperbolicRadius, double &euclideanCenter, double &euclideanRadius);
	static double hyperbolicRadiusToEuclidean(double hyperbolicRadius);
	static double EuclideanRadiusToHyperbolic(double EuclideanRadius);
	static double hyperbolicSpaceInEuclideanCircle(double r_c, double d_c, double R);
};
}

#endif /* HYPERBOLICSPACE_H_ */
