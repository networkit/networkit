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
	HyperbolicSpace(double R);
	static void fillPoints(vector<double> * angles, vector<double> * radii, double stretch, double alpha);
	static double getHyperbolicDistance(double firstangle, double firstR, double secondangle, double secondR);
	static double getHyperbolicDistance(Point2D<double> a, Point2D<double> b);
	double getRadius();
	static Point2D<double> polarToCartesian(double phi, double r);
	static void cartesianToPolar(Point2D<double> a, double &phi, double &r);
	static void getEuclideanCircle(Point2D<double> hyperbolicCenter, double hyperbolicRadius, Point2D<double> &euclideanCenter, double &euclideanRadius);
	static void getEuclideanCircle(double r_h, double hyperbolicRadius, double &euclideanCenter, double &euclideanRadius);
	static double hyperbolicRadiusToEuclidean(double hyperbolicRadius);
	static double EuclideanRadiusToHyperbolic(double EuclideanRadius);
	static double hyperbolicSpaceInEuclideanCircle(double r_c, double d_c, double R);


private:
	double radius;
	//one could add some caching here, but this should be done properly in an object. Well, or maybe not.
	static double lastR;
	static double coshlastR;
	static double sinhlastR;
};
}

#endif /* HYPERBOLICSPACE_H_ */
