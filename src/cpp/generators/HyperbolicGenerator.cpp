/*
 * HyperbolicGenerator.cpp
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#include <cstdlib>
#include <random>
#include <math.h>
#include <assert.h>

#include "../graph/GraphBuilder.h"
#include "HyperbolicGenerator.h"
#include "Quadtree/Quadtree.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/ProgressMeter.h"

namespace NetworKit {


HyperbolicGenerator::HyperbolicGenerator() {
	stretch = 1;
	alpha = 1;
	factor = 1;
	nodeCount = 10000;
}

HyperbolicGenerator::HyperbolicGenerator(count n, double distanceFactor, double alpha, double stretchradius) {
	nodeCount = n;
	stretch = stretchradius;
	factor = distanceFactor;
	this->alpha = alpha;
}

HyperbolicGenerator::~HyperbolicGenerator() {
	// TODO Auto-generated destructor stub
}

Graph HyperbolicGenerator::generate() {
	return generate(nodeCount, factor, alpha, stretch);
}

Graph HyperbolicGenerator::generate(count n, double distanceFactor, double alpha, double stretchradius) {
	double R = stretchradius*acosh((double)n/(2*M_PI)+1);
	vector<double> angles(n);
	vector<double> radii(n);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	HyperbolicSpace::fillPoints(&angles, &radii, stretchradius, alpha);
	INFO("Generated Points");
	return generate(&angles, &radii, r, R*distanceFactor);
}

double HyperbolicGenerator::expectedNumberOfEdges(count n, double distanceFactor) {
	double R = acosh((double)n/(2*M_PI)+1);
	double euR = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	double radius = R*distanceFactor;
	double epsilon = 0.00001;
	double upperlimit =-log(1-euR);

	//numeric integration
	double stepsize = 0.01;
	double result;
	for (double step = 0; step < upperlimit; step += stepsize) {
		double r_a = 1-(exp(-step));
		double r_a_next = 1-(exp(-(step+stepsize)));
		double stSq = r_a*r_a;
		double r_c = (2*r_a)/((1-stSq)*(cosh(radius)-1+2/(1-stSq)));
		double euRadius = sqrt(r_c*r_c - (2*stSq - (1-stSq)*(cosh(radius)-1))/((1-stSq)*(cosh(radius)-1+2/(1-stSq))));
		double circleSpace;
		if (r_c*r_c-euRadius*euRadius+1 <= 2*r_c) {
			circleSpace = 0;
		} else {
			circleSpace = HyperbolicSpace::hyperbolicSpaceInEuclideanCircle(r_c,euRadius,euR);
		}
		if (!(circleSpace == circleSpace)) circleSpace = 0;//careful, this will behave differently in Debug mode and optimized
		assert(circleSpace <= n+epsilon);
		double outerSpace = 2*M_PI*(cosh(HyperbolicSpace::EuclideanRadiusToHyperbolic(std::min(r_a_next,euR)))-1);
		double innerSpace = 2*M_PI*(cosh(HyperbolicSpace::EuclideanRadiusToHyperbolic(r_a))-1);
		result += (outerSpace-innerSpace)*circleSpace;
		assert(result == result);
	}
	return result / 2;
}

Graph HyperbolicGenerator::generate(vector<double> * angles, vector<double> * radii, double R, double thresholdDistance) {
	index n = angles->size();
	assert(radii->size() == n);
	Quadtree<index> quad(R);
	GraphBuilder result(n, false, false, true);
	for (index i = 0; i < n; i++) {
		assert(radii->at(i) < R);
		quad.addContent(i, angles->at(i), radii->at(i));
	}
	INFO("Filled Quadtree");

	Aux::ProgressMeter progress(n, 1000);
	#pragma omp parallel for schedule(dynamic, 1000)
	for (index i = 0; i < n; i++) {
			vector<index> near = quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles->at(i), radii->at(i)), thresholdDistance);
			for (index j : near) {
				if (i != j) {
						result.addEdge(i,j);
				}
			}

			if (i % 1000 == 0) {
				#pragma omp critical (progress)//that doesn't make any sense, creating the block every time and only printing every 200th iterations
				{
					progress.signal(i);
				}
			}
		}

	return result.toGraph(true);
}


}
