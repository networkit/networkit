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
	GraphGenerator graphGen;
	Graph G = graphGen.makeRandomGraph(n, 0.2);
	G.initCoordinates();

	// create clustering
	ClusteringGenerator clusteringGen;
	Clustering zeta = clusteringGen.makeRandomClustering(G, 3);

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
	count n = 50;
	GraphGenerator graphGen;
	Graph G = graphGen.makeRandomGraph(n, 0.2);
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
	Clustering zeta = clusteringGen.makeRandomClustering(G, 3);

	// draw (independent of clustering)
	std::vector<float> bottomLeft = {0.0, 0.0};
	std::vector<float> topRight = {1000.0, 1000.0};
	Point<float> bl(bottomLeft);
	Point<float> tr(topRight);
	ForceDirected layouter(bl, tr);
	layouter.draw(G);

	// write graph to file
	PostscriptWriter psWriter(G);
	psWriter.write(zeta, "testGraph.eps");
}




} /* namespace NetworKit */

