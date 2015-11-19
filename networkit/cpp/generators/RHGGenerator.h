#ifndef RHGGenerator_H_
#define RHGGenerator_H_

#include <vector>
#include "../geometric/HyperbolicSpace.h"
#include "StaticGraphGenerator.h"
#include "../auxiliary/Timer.h"
#include <assert.h>
#include <cmath>
#include <tuple>
#include <utility>
#include <iostream>


#include "../auxiliary/Log.h"

using std::abs;
using std::max;
using std::vector;

namespace NetworKit {

	typedef uint64_t index; // more expressive name for an index into an array
	typedef uint64_t count; // more expressive name for an integer quantity
	typedef index node; // node indices are 0-based

	class RHGGenerator: public NetworKit::StaticGraphGenerator {
	public:
		RHGGenerator();

		/**
		* @param[in] n Number of nodes
		*/
		RHGGenerator(count n);

		/**
		* @param[in] n Number of nodes
		* @param[in] m Target number of edges
		*/
		RHGGenerator(count n, double avgDegree=6, double exp=3);


		Graph generate(const vector<double> &angles, const vector<double> &radii, const vector<vector<Point2D<double>>> &bands, const vector<double> &bandRadius, double thresholdDistance);

		/**
		* @param[in] angles Pointer to angles of node positions
		* @param[in] radii Pointer to radii of node positions
		* @param[in] r radius of poincare disk to place nodes in
		* @param[in] thresholdDistance Edges are added for nodes closer to each other than this threshold
		* @return Graph to be generated according to parameters
		*/
		Graph generate(const vector<double> &angles, const vector<double> &radii, double thresholdDistance);

		/**
		* @param[in] n Number of nodes
		* @param[in] factor Size of neighborhood radius. If factor=1, radius = R
		* @param[in] alpha Dispersion parameter, default=1
		* @param[in] stretchradius Stretching the hyperbolic disk results in thinner graphs, default=1
		* @return Graph to be generated according to parameters
		*/
		Graph generate(count n, double distanceFactor=1, double alpha=1, double stretchradius = 1);

		/**
		* @return Graph to be generated according to parameters specified in constructor.
		*/
		Graph generate();

		vector<double> getElapsedMilliseconds() {
			vector<double> result(threadtimers.size());
			for (index i = 0; i < result.size(); i++) {
				result[i] = threadtimers[i].elapsedMilliseconds();
			}
			return result;
		}

	private:
		/**
		* graph parameters
		*/
		count nodeCount;
		double stretch;
		double factor;
		double alpha;

		static const bool directSwap = false;

		/**
		* times
		*/
		vector<Aux::Timer> threadtimers;

		/**
		* helper methods
		*/

		void getBandRadius(int n, vector<double> &bandRadius, double thresholdDistance, double seriesRatio = 0.5) {
			/*
			* We asumme band differences form a geometric series.
			* Thus, there is a constant ratio(r) between band length differences
			* i.e c2-c1/c1-c0 = c3-c2/c2-c1 = r
			*/
			bandRadius.push_back(0);
			double a = thresholdDistance*(1-seriesRatio)/(1-pow(seriesRatio, log(n)));

			for (int i = 1; i < ceil(log(n)); i++){
				double c_i = a*((1-pow(seriesRatio, i))/(1-seriesRatio));
				bandRadius.push_back(c_i);
			}
			bandRadius.push_back(thresholdDistance);
		}

		std::tuple<double, double> getMinMaxTheta(double angle, double radius, double cLow, double thresholdDistance){
			/*
			Calculates the angles that are enclosing the intersection of the
			hyperbolic disk that is around point v and the bands.
			Calculation is as follows:
			1. For the most inner band, return [0, 2pi]
			2. For other bands, consider the point P which lies on the tangent from origin to the disk of point v.
			Its radial coordinates would be(cHigh, point[1]+deltaTheta). We're looking for the deltaTheta.
			We know the distance from point v to P is R. Thus, we can solve the hyperbolic distance of (v, P)
			for deltaTheta. Then, thetaMax is simply point[1] + deltaTheta and thetaMin is point[1] - deltaTheta
			*/

			//Most innerband is defined by cLow = 0
			double minTheta, maxTheta;
			if (cLow == 0)
			return std::make_tuple(0, 2* M_PI);

			double a = acos((cosh(radius)*cosh(cLow) - cosh(thresholdDistance))/(sinh(radius)*sinh(cLow)));
			maxTheta = angle + a;
			minTheta = angle - a;
			return std::make_tuple(minTheta, maxTheta);
		}


		void getPointsWithinAngles(double minTheta, double maxTheta, const vector<Point2D<double>> &band, vector<double> &bandAngles, vector<Point2D<double>> &slab){
			/**
			Returns the list of points, w, that lies within minTheta and maxTheta
			in the supplied band(That area is called as slab)
			*/
			//TODO: There should be a better way to write the whole thing. Find it.

			std::vector<double>::iterator low;
			std::vector<double>::iterator high;

			//Case 1: We do not have overlap 2pi, simply put all the points between min and max to the list
			if(maxTheta <= 2*M_PI && minTheta >= 0){
				low = std::lower_bound(bandAngles.begin(), bandAngles.end(), minTheta);
				high = std::upper_bound(bandAngles.begin(), bandAngles.end(), maxTheta);
				std::vector<Point2D<double>>::const_iterator first = band.begin() + (low - bandAngles.begin());
				std::vector<Point2D<double>>::const_iterator last = band.begin() + (high - bandAngles.begin());
				//Q: Does this operation increases the complexity ? It is linear in times of high - low
				slab.insert(slab.end(), first, last);
				return;
			}
			//Case 2: We have 'forward' overlap at 2pi, that is maxTheta > 2pi
			else if (maxTheta > 2*M_PI){
				//1. Get points from minTheta to 2pi
				low = std::lower_bound(bandAngles.begin(), bandAngles.end(), minTheta);
				high = std::upper_bound(bandAngles.begin(), bandAngles.end(), 2*M_PI);
				std::vector<Point2D<double>>::const_iterator first = band.begin() + (low - bandAngles.begin());
				std::vector<Point2D<double>>::const_iterator last = band.begin() + (high - bandAngles.begin());
				slab.insert(slab.end(), first, last);

				//2. Get points from 0 to maxTheta%2pi
				low = std::lower_bound(bandAngles.begin(), bandAngles.end(), 0);
				maxTheta = fmod(maxTheta, (2*M_PI));
				high = std::upper_bound(bandAngles.begin(), bandAngles.end(), maxTheta);
				std::vector<Point2D<double>>::const_iterator first2 = band.begin() + (low - bandAngles.begin());
				std::vector<Point2D<double>>::const_iterator last2 = band.begin() + (high - bandAngles.begin());
				slab.insert(slab.end(), first2, last2);

				return;
			}
			//Case 3: We have 'backward' overlap at 2pi, that is minTheta < 0
			else if (minTheta < 0){
				//1. Get points from 2pi - minTheta to 2pi
				minTheta = (2*M_PI) - minTheta;
				low = std::lower_bound(bandAngles.begin(), bandAngles.end(), minTheta);
				high = std::upper_bound(bandAngles.begin(), bandAngles.end(), 2*M_PI);
				std::vector<Point2D<double>>::const_iterator first = band.begin() + (low - bandAngles.begin());
				std::vector<Point2D<double>>::const_iterator last = band.begin() + (high - bandAngles.begin());
				slab.insert(slab.end(), first, last);
				//2. Get points from 0 to maxTheta
				low = std::lower_bound(bandAngles.begin(), bandAngles.end(), 0);
				high = std::upper_bound(bandAngles.begin(), bandAngles.end(), maxTheta);
				std::vector<Point2D<double>>::const_iterator first2 = band.begin() + (low - bandAngles.begin());
				std::vector<Point2D<double>>::const_iterator last2 = band.begin() + (high - bandAngles.begin());
				slab.insert(slab.end(), first2, last2);
				return;
			}
		}


		double getHyperbolicDistance(Point2D<double> u, Point2D<double> v) {
			/* Returns the hyperbolic distance between points u and v
			* 2010 paper, eqn: 5
			*/
			double deltaTheta = M_PI - abs(M_PI-abs(u.getX() - v.getX()));
			double distance = acosh(cosh(u.getY())*cosh(v.getY()) - sinh(u.getY())*sinh(v.getY())*cos(deltaTheta));
			return distance;
		}


		/*
		* Modified methods from HyperbolicSpace.h
		* Native representation of the hyperbolic disk
		* is used instead of the Poincare
		*
		*/
		void fillPoints(vector<double> &angles, vector<double> &radii, double stretch, double alpha) {
			uint64_t n = radii.size();
			double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(n);
			assert(angles.size() == n);
			double minPhi = 0;
			double maxPhi = 2*M_PI;
			double minR = 0;
			double maxR = R;

			double mincdf = cosh(alpha*minR);
			double maxcdf = cosh(alpha*maxR);
			std::uniform_real_distribution<double> phidist{minPhi, maxPhi};
			std::uniform_real_distribution<double> rdist{mincdf, maxcdf};
			//do not transform to Poincare
			double r = maxR;

			for (uint64_t i = 0; i < n; i++) {
				angles[i] = phidist(Aux::Random::getURNG());
				/**
				* for the radial coordinate distribution, I took the probability density from Greedy Forwarding in Dynamic Scale-Free Networks Embedded in Hyperbolic Metric Spaces
				* f (r) = sinh r/(cosh R − 1)
				* \int sinh = cosh+const
				*/

				double random = rdist(Aux::Random::getURNG());
				double radius = (acosh(random)/alpha);
				//assert(radius < maxR);
				radii[i] = radius;
				assert(radii[i] <= r);
				if (radii[i] == r) radii[i] = std::nextafter(radii[i], 0);
				assert(radii[i] < r);
			}
		}
	};
}
#endif /* RHGGenerator_H_ */
