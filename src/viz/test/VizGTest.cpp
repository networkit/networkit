/*
 * PostscriptWriterGTest.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#include "VizGTest.h"


namespace NetworKit {

VizGTest::VizGTest() {

}

VizGTest::~VizGTest() {

}


TEST_F(VizGTest, testPostscriptWriter) {
	// create graph
	count n = 60;
	count numClusters = 3;
	double pin = 0.35;
	double pout = 0.05;

	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, numClusters, pin, pout);
	G.initCoordinates();

	// create clustering
	ClusteringGenerator clusteringGen;
	Clustering zeta = clusteringGen.makeRandomClustering(G, numClusters);

	// create coordinates
	G.forNodes([&](node u) {
		float x = (float) drand48();
		float y = (float) drand48();
		G.setCoordinate(u, 0, x);
		G.setCoordinate(u, 1, y);
	});

	// write graph to file
	PostscriptWriter psWriter(G);
	psWriter.write(zeta, "testGraph.eps");
}


TEST_F(VizGTest, testLayouter) {
	// create graph
	count n = 50;
	count numClusters = 2;
	double pin = 0.35;
	double pout = 0.05;

	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, numClusters, pin, pout);
	G.initCoordinates();
	INFO("Number of edges: " << G.numberOfEdges());

	// create clustering
	ClusteringGenerator clusteringGen;

	// write graph to file
	PostscriptWriter psWriter(G, true);
	psWriter.write("output/testGraph.eps");

	// draw (independent of clustering) and write again
	std::vector<float> bottomLeft = {0.0, 0.0};
	std::vector<float> topRight = {1.0, 1.0};
	Point<float> bl(bottomLeft);
	Point<float> tr(topRight);

	ForceDirected fdLayouter(bl, tr);
	fdLayouter.draw(G);
	PostscriptWriter psWriter2(G, true);
	psWriter2.write("output/testForceGraph.eps");

	MaxentStress msLayouter(bl, tr);
	msLayouter.draw(G);
	PostscriptWriter psWriter3(G, true);
	psWriter3.write("output/testMaxentGraph.eps");
}



} /* namespace NetworKit */


#endif /*NOGTEST */


