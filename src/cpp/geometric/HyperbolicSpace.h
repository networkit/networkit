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
#include "../viz/Point.h"

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
	static double getHyperbolicDistance(Point<double> a, Point<double> b);
	double getRadius();
	static double cross(Point<double> a, Point<double> b);
	static Point<double> intersect(Point<double> q, Point<double> s, Point<double> p, Point<double> r);
	static Point<double> circleCenter(Point<double> a, Point<double> b, Point<double> c);
	static Point<double> mirrorOnCircle(Point<double> a, Point<double> circleCenter, double radius);
	static bool isBelowArc(Point<double> query, Point<double> a, Point<double> b, double radius);
	static Point<double> polarToCartesian(double phi, double r);
	static void cartesianToPolar(Point<double> a, double &phi, double &r);
	static Point<double> orth(Point<double> a);
	static void getTransmutationCircle(Point<double> source, Point<double> target, double minRadius, Point<double> &circleCenter, double &circleRadius);
	static double hyperbolicDistanceToArc(Point<double> query, Point<double> a, Point<double> b, double R);
	static Point<double> getPointOnHyperbolicCircle(Point<double> hyperbolicCenter, double radius);
	static void getEuclideanCircle(Point<double> hyperbolicCenter, Point<double> hyperbolicRadius, Point<double> &euclideanCenter, double &euclideanRadius);
	static double hyperbolicRadiusToEuclidean(double hyperbolicRadius);
	static double EuclideanRadiusToHyperbolic(double EuclideanRadius);


private:
	double radius;
	//one could add some caching here, but this should be done properly in an object. Well, or maybe not.
	static double lastR;
	static double coshlastR;
	static double sinhlastR;
};
}

#endif /* HYPERBOLICSPACE_H_ */
