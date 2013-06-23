/*
 * DistMeasureTest.cpp
 *
 *  Created on: Jun 11, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#include "DistMeasuresGTest.h"

namespace NetworKit {

DistMeasuresGTest::DistMeasuresGTest() {

}

DistMeasuresGTest::~DistMeasuresGTest() {

}

TEST_F(DistMeasuresGTest, testAlgebraicDistances) {
	GraphGenerator gen;
	Graph G = gen.makeCompleteGraph(42);

	count numSystems = 2;
	count numIterations = 200;
	double omega = 0.5;
	NodeDistance* ad = new AlgebraicDistance(G, numSystems, numIterations, omega);
	ad->preprocess();

	double adSum = 0.0;
	G.forNodePairs([&](node u, node v){
		adSum += ad->distance(u, v);
	});

	INFO("sum of algebraic distances: " << adSum);
	EXPECT_GE(1e-12, adSum) << "algebraic distances should converge towards zero";

}

TEST_F(DistMeasuresGTest, tryAlgebraicDistancesVisualization) {
#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64

	DibapGraphReader reader;
	std::string path = "input/airfoil1.gi";
	Graph g = reader.read(path);

	// init algebraic distances and preprocess
//	count n = g.numberOfNodes();
	count numSystems = 8;
	count numIterations = 25;
	double omega = 0.5;
//	const count norm = 2;
	AlgebraicDistances ad(g);
	ad.preprocess(numSystems, numIterations, omega);

	PostscriptWriter psWriter(g);
	psWriter.writeAlgebraicDistances("output/airfoil1-ad-test.eps", ad);

#if 0
	// *** produce 2-clustering according to AD ***/

	// determine two vertices with highest distance
	index best0 = 0;
	index best1 = 1;
	double maxDist = ad.algdist(best0, best1, norm);
	g.forNodePairs([&](node i, node j) {
		double dist = ad.algdist(i, j , norm);
		if (dist > maxDist) {
			best0 = i;
			best1 = j;
			maxDist = dist;
		}
	});
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
#endif

#endif // non-Windows test
}


} /* namespace NetworKit */

#endif /*NOGTEST */
