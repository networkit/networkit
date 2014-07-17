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

HyperbolicGenerator::HyperbolicGenerator(count n, double stretchradius, double distanceFactor) {
	nodeCount = n;
	stretch = stretchradius;
	factor = distanceFactor;
}

HyperbolicGenerator::~HyperbolicGenerator() {
	// TODO Auto-generated destructor stub
}

Graph HyperbolicGenerator::generate() {
	return generate(nodeCount, stretch, factor);
}

Graph HyperbolicGenerator::generate(count n, double stretchradius, double distanceFactor) {
	double R = stretchradius*acosh((double)n/(2*M_PI)+1);
	vector<double> angles(n);
	vector<double> radii(n);
	double rad_nom = (cosh(R)-1);
	double rad_denom = (cosh(R)+1);
	double r = sqrt(rad_nom/rad_denom);
	HyperbolicSpace::fillPoints(&angles, &radii, stretchradius, 1);
	INFO("Generated Points");
	return generate(&angles, &radii, r, R*distanceFactor);
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
	INFO("Filled Quadtree");

	Aux::ProgressMeter progress(n, 200);
	#pragma omp parallel for
	for (index i = 0; i < n; i++) {
			vector<index> near = quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles->at(i), radii->at(i)), thresholdDistance);
			TRACE("gathered points for ", i, " adding edges");
			for (index j : near) {
				if (i != j) {//we only want to add the edges once for each pair
						result.addHalfEdge(i,j);
				}
			}
			TRACE("added edges");

			if (i % 200 == 0) {
				TRACE("giving progress signal");
				#pragma omp critical (progress)//that doesn't make any sense, creating the block every time and only printing every 200th iterations
				{
					progress.signal(i);
					TRACE("total edges:", result.totalEdgeWeight());
				}
			}
		}

	return result;
}


}
