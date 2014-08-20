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
