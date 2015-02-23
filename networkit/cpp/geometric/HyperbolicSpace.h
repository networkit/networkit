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
#include "../auxiliary/Random.h"

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
	static void fillPoints(vector<double> &angles, vector<double> &radii, double stretch, double alpha);

	/**
	 * @param angles empty vector to hold angular coordinates of generated points
	 * @param radii empty vector to hold radial coordinates of generated points
	 * @param stretch multiplier for the radius of the hyperbolic disk
	 * @param alpha dispersion parameter for the node positions
	 *
	 * Fill preallocated vectors with randomly sampled points in the poincare disk
	 */
	static void fillPoints(vector<double> &angles, vector<double> &radii, double minPhi, double maxPhi, double minR, double maxR, double alpha);

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
	 * @param EuclideanRadius Radial coordinate of a point in the Poincare model
	 * @return radial coordinate in the hyperbolic plane
	 */
	static double EuclideanRadiusToHyperbolic(double EuclideanRadius);

	/**
	 * @param area The area of the hyperbolic circle
	 * @param return Radius of a hyperbolic circle with the given area
	 */
	static inline double hyperbolicAreaToRadius(double area) {
		return acosh(area/(2*M_PI)+1);
	}

	static inline double radiusToHyperbolicArea(double radius) {
			return  2*M_PI*(cosh(radius)-1);
	}

	static double getTargetRadius(double n, double m, double alpha, double T) {
		double result;
		double plexp = 2*alpha+1;
		double avgDegree = m/n;
		double gammaratio = ((plexp-2)/(plexp-1));
		if (T == 0) {
			result = 2*log(8*n / (M_PI*(m/n)*2));
		} else {
			double beta = 1/T;
			if (T < 1){//cold regime
				double Iinv = ((beta/M_PI)*sin(M_PI/beta));
				double v = (avgDegree*Iinv)*(M_PI/2)*gammaratio*gammaratio;
				result = 2*log(n / v);
			} else {//hot regime
				double v = avgDegree*(1-beta)*pow((M_PI/2), beta)*gammaratio*gammaratio;
				result = 2*log(n/v)/beta;
			}
		}
		return result;
	}

	static inline double effectiveAreaInCell(double minPhi, double maxPhi, double minR, double maxR, double alpha) {
		double deltaPhi = maxPhi - minPhi;
		assert(deltaPhi >= 0);
		assert(deltaPhi <= 2*M_PI);
		double ringArea = radiusToHyperbolicArea(alpha*EuclideanRadiusToHyperbolic(maxR)) - radiusToHyperbolicArea(alpha*EuclideanRadiusToHyperbolic(minR));
		return ringArea*(deltaPhi/(2*M_PI));
	}

	/**
	 * @param r_c radial coordinate of the circle center
	 * @param d_c radius of the Euclidean circle
	 * @param R radius of the hyperbolic base disk
	 */
	static double hyperbolicSpaceInEuclideanCircle(double r_c, double d_c, double R);

	/**
	 *
	 */
	static double maxRinSlice(double minPhi, double maxPhi, double phi_c, double r_c, double euRadius);
};
}

#endif /* HYPERBOLICSPACE_H_ */
