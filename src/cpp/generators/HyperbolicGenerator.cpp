/*
 * HyperbolicGenerator.cpp
 *
 *  Created on: 20.05.2014
 *      Author: moritz
 */

#include <cstdlib>
#include <random>
#include <math.h>
#include <assert.h>

#include "HyperbolicGenerator.h"
#include "Quadtree/Quadtree.h"
#include "../auxiliary/Random.h"

namespace NetworKit {


HyperbolicGenerator::HyperbolicGenerator() {
	stretch = 1;
	nodeCount = 10000;
}

HyperbolicGenerator::HyperbolicGenerator(count n, double stretchradius) {
	nodeCount = n;
	stretch = stretchradius;
}

HyperbolicGenerator::~HyperbolicGenerator() {
	// TODO Auto-generated destructor stub
}

Graph HyperbolicGenerator::generate() {
	return generate(nodeCount, stretch);
}

Graph HyperbolicGenerator::generate(count n, double stretchradius) {
	double R = stretchradius*acosh((double)n/(2*M_PI)+1);
	HyperbolicSpace space(R);
	Graph result(n, false);
	vector<double> angles(n);
	vector<double> radii(n);
	Quadtree<index> quad(R);

	for (index i = 0; i < n; i++) {
		angles[i] = ((double)  rand() / (double) RAND_MAX)*2*M_PI;
		/**
		 * for the radial coordinate distribution, I took the probability density from Greedy Forwarding in Dynamic Scale-Free Networks Embedded in Hyperbolic Metric Spaces
		 * f (r) = sinh r/(cosh R − 1)
		 * \int sinh = cosh+const
		 */

		//double denominator = cosh(R)-1;//not needed because it is a common factor
		double maxcdf = cosh(R);
		double random = Aux::Random::real(1, maxcdf);
		radii[i] = acosh(random);
		assert(radii[i] <= R);
		quad.addContent(i, angles[i], radii[i]);
	}

	for (index i = 0; i < n; i++) {
			vector<index> near = quad.getCloseElements(angles[i], radii[i], R);
			for (index j : near) {
				if (i < j) {//we only want to add the edges once for each pair
					result.addEdge(i,j);
				}
			}
		}
	return result;
}
}
