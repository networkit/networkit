/*
 * DistMeasureTest.cpp
 *
 *  Created on: Jun 11, 2013
 *      Author: Henning
 */

#include "DistMeasureTest.h"

namespace NetworKit {

DistMeasureTest::DistMeasureTest() {

}

DistMeasureTest::~DistMeasureTest() {

}

TEST_F(DistMeasureTest, tryAlgebraicDistances) {
	METISGraphReader reader;
	Graph g = reader.read("input/airfoil1.graph");
	std::cout.precision(5);

	// init algebraic distances and preprocess
	count n = g.numberOfNodes();
	count numSystems = 5;
	count numIterations = 20;
	double omega = 0.5;
	const count norm = 2;
	AlgebraicDistances ad(g);
	ad.preprocess(numSystems, numIterations, omega);

#if 0
	// print pairwise distances for a subset of vertices
	for (index i = 0; i < 20; ++i) {
		for (index j = 0; j < i; ++j) {
			printf("%.4f  ", ad.algdist(i, j, norm));
		}
		std::cout << std::endl;
	}
#endif

	// write color-coded clustering according to AD

	// determine two vertices with highest distance
	index best0 = 0;
	index best1 = 0;
	double maxDist = ad.algdist(best0, best1, norm);

	for (index i = 0; i < n; ++i) {
		for (index j = 0; j < i; ++j) {
			double dist = ad.algdist(i, j , norm);
			if (dist > maxDist) {
				best0 = i;
				best1 = j;
				maxDist = dist;
			}
		}
	}
	DEBUG("Finished with maxDist computation");

	// assign remaining vertices to one of the two according to the smaller distance
	Clustering clustering(n);
	clustering[best0] = 0;
	clustering[best1] = 1;
	g.forNodes([&](node u) {
		clustering[u] = (ad.algdist(u, best0) < ad.algdist(u, best1)) ? (0) : (1);
	});

	Modularity modScore;

	DEBUG("Finished with assignment to clusters, mod score: " << modScore.getQuality(clustering, g));


//	PostscriptWriter psWriter(g);
//	psWriter.write(clustering, "output/airfoil1-ad-test.eps");
}


} /* namespace NetworKit */
