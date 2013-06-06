/*
 * PostscriptWriterGTest.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#include "VizGTest.h"


namespace NetworKit {

VizGTest::VizGTest() {
	// TODO Auto-generated constructor stub

}

VizGTest::~VizGTest() {
	// TODO Auto-generated destructor stub
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


TEST_F(VizGTest, testForceDirectedLayouter) {
	// create graph
	count n = 30;
	count numClusters = 2;
	double pin = 0.35;
	double pout = 0.05;

	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, numClusters, pin, pout);
	G.initCoordinates();

	// create coordinates
	G.forNodes([&](node u) {
		float x = (float) drand48();
		float y = (float) drand48();
		G.setCoordinate(u, 0, x);
		G.setCoordinate(u, 1, y);
	});

	// create clustering
	ClusteringGenerator clusteringGen;
	Clustering zeta = clusteringGen.makeOneClustering(G);

	// write graph to file
	PostscriptWriter psWriter(G);
	psWriter.write(zeta, "testGraph.eps");

	// draw (independent of clustering) and write again
	std::vector<float> bottomLeft = {0.0, 0.0};
	std::vector<float> topRight = {1.0, 1.0};
	Point<float> bl(bottomLeft);
	Point<float> tr(topRight);
	ForceDirected layouter(bl, tr);
	layouter.draw(G);
	psWriter.write(zeta, "testForceGraph.eps");
}




} /* namespace NetworKit */

