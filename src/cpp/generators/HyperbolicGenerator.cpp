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

#include "HyperbolicGenerator.h"
#include "Quadtree/Quadtree.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/ProgressMeter.h"

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
	vector<double> angles(n);
	vector<double> radii(n);

	HyperbolicSpace::fillPoints(&angles, &radii, stretchradius, 1);
	INFO("Generated Points");
	return generate(&angles, &radii, R, R);
}

Graph HyperbolicGenerator::generate(vector<double> * angles, vector<double> * radii, double R, double thresholdDistance) {
	index n = angles->size();
	assert(radii->size() == n);
	Quadtree<index> quad(R);
	Graph result(n, false);
	for (index i = 0; i < n; i++) {
		assert(radii->at(i) < R);
		quad.addContent(i, angles->at(i), radii->at(i));
	}
	Aux::ProgressMeter progress(n, 500);
	INFO("Filled Quadtree");
	#pragma omp parallel for
	for (index i = 0; i < n; i++) {
			vector<index> near = quad.getCloseElements(angles->at(i), radii->at(i), R);
			for (index j : near) {
				if (i < j) {//we only want to add the edges once for each pair
					#pragma omp critical
					{
					result.addEdge(i,j);
					}
				}
			}
			#pragma omp critical
			progress.signal(i);
		}

	return result;
}


}
