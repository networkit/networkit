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
	count n = 150;
	GraphGenerator graphGen;
	Graph G = graphGen.makeRandomGraph(n, 0.2);

	// create clustering
	ClusteringGenerator clusteringGen;
	Clustering zeta = clusteringGen.makeRandomClustering(G, 3);

	// create coordinates
	G.coordinates.init(n, 2);
	G.forNodes([&](node u) {
		float x = (float) drand48();
		float y = (float) drand48();
		G.coordinates.setCoordinate(u, 0, x);
		G.coordinates.setCoordinate(u, 1, y);
	});

	// write graph to file
	PostscriptWriter psWriter(G);
	psWriter.write(zeta, "testGraph.eps");
}



} /* namespace NetworKit */

