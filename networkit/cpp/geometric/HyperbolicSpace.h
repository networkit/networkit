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
	 * @return distance between two points in the poincare metric
	 */
	static double poincareMetric(double firstangle, double firstR, double secondangle, double secondR);

	/**
	 * @param a first point in cartesian coordinates
	 * @param b second point in cartesian coordinates
	 *
	 * @return distance between a and b in the poincare metric
	 */
	static double poincareMetric(Point2D<double> a, Point2D<double> b);

	/**
	 * @param phi angular coordinate of point
	 * @param r radial coordinate of point
	 *
	 * @return cartesian coordinates represented by phi and r
	 */
	static Point2D<double> polarToCartesian(double phi, double r);

	/**
	 * Convenience function for visualizations which expect coordinates as map<index,Point<float> >
	 */
	static std::map<index, Point<float> > polarToCartesian(vector<double> &angles, vector<double> &radii);

	/**
	 * @param a cartesian coordinates
	 * @param phi empty double value to receive angular coordinate
	 * @param r empty double value to receive radial coordinate
	 */
	static void cartesianToPolar(Point2D<double> a, double &phi, double &r);

	/**
	 * Converts a hyperbolic circle to a Euclidean circle
	 *
	 * @param hyperbolicCenter center of the hyperbolic circle, given in cartesian coordinates within the poincare disk
	 * @param hyperbolicRadius radius of the hyperbolic circle
	 * @param euclideanCenter point to receive the center of the Euclidean circle, given in cartesian coordinates
	 * @param euclidenRadius double to receive the radius of the Euclidean circle
	 */
	static void getEuclideanCircle(Point2D<double> hyperbolicCenter, double hyperbolicRadius, Point2D<double> &euclideanCenter, double &euclideanRadius);
	static void getEuclideanCircle(double r_h, double hyperbolicRadius, double &radialCoordOfEuclideanCenter, double &euclideanRadius);

	/**
	 * Project radial coordinates of the hyperbolic plane into the Poincare disk model
	 *
	 * @param hyperbolicRadius radial coordinate of a point in the native hyperbolic disc
	 * @return radial coordinate in the Poincare model
	 */
	static double hyperbolicRadiusToEuclidean(double hyperbolicRadius);

	/**
	 * Project radial coordinates of the Poincare model into the hyperbolic plane
	 *
	 * @param EuclideanRadius radial coordinate of a point in the Poincare model
	 * @return radial coordinate in the hyperbolic plane
	 */
	static double EuclideanRadiusToHyperbolic(double EuclideanRadius);

	/**
	 * @param r_c radial coordinate of the circle center
	 * @param d_c radius of the Euclidean circle
	 * @param R radius of the hyperbolic base disk
	 */
	static double hyperbolicSpaceInEuclideanCircle(double r_c, double d_c, double R);
};
}

#endif /* HYPERBOLICSPACE_H_ */
