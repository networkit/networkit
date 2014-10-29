/*
 * HyperbolicSpace.h
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef HYPERBOLICSPACE_H_
#define HYPERBOLICSPACE_H_

#include <vector>
#include <map>
#include <cmath>
#include "Point2D.h"
#include "../viz/Point.h"

using std::vector;
using std::abs;

namespace NetworKit {

class HyperbolicSpace {
public:
	HyperbolicSpace();
	virtual ~HyperbolicSpace();
	/**
	 * @param angles empty vector to hold angular coordinates of generated points
	 * @param radii empty vector to hold radial coordinates of generated points
	 * @param stretch multiplier for the radius of the hyperbolic disk
	 * @param alpha dispersion parameter for the node positions
	 *
	 * Fill preallocated vectors with randomly sampled points in the poincare disk
	 */
	static void fillPoints(vector<double> * angles, vector<double> * radii, double stretch, double alpha);

	/**
	 * @param firstangle angular coordinate of the first point
	 * @param firstR radial coordiante of the first point
	 * @param secondangle angular coordinate of the second point
	 * @param secondR radial coordinate of the second point
	 *
	 * @return distance in the poincare metric
	 */
	static double poincareMetric(double firstangle, double firstR, double secondangle, double secondR);
	static double poincareMetric(Point2D<double> a, Point2D<double> b);
	static Point2D<double> polarToCartesian(double phi, double r);
	/**
	 * Convenience function for visualizations which expect coordinates as map<index,Point<float> >
	 */
	static std::map<index, Point<float> > polarToCartesian(vector<double> &angles, vector<double> &radii);
	static void cartesianToPolar(Point2D<double> a, double &phi, double &r);
	static void getEuclideanCircle(Point2D<double> hyperbolicCenter, double hyperbolicRadius, Point2D<double> &euclideanCenter, double &euclideanRadius);
	static void getEuclideanCircle(double r_h, double hyperbolicRadius, double &euclideanCenter, double &euclideanRadius);
	static double hyperbolicRadiusToEuclidean(double hyperbolicRadius);
	static double EuclideanRadiusToHyperbolic(double EuclideanRadius);
	static double hyperbolicSpaceInEuclideanCircle(double r_c, double d_c, double R);
};
}

#endif /* HYPERBOLICSPACE_H_ */
