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
		Point<float> p(drand48(), drand48());
		G.setCoordinate(u, p);
	});

	// write graph to file
	PostscriptWriter psWriter(G);
	psWriter.write(zeta, "testGraph.eps");
}

TEST_F(VizGTest, testFRLayouter) {
 	// create graph
 	count n = 120;
 	count numClusters = 3;
 	double pin = 0.175;
 	double pout = 0.005;

 	GraphGenerator graphGen;
 	Graph G = graphGen.makeClusteredRandomGraph(n, numClusters, pin, pout);
 	G.initCoordinates();
 	INFO("Number of edges: ", G.numberOfEdges());

 	// draw (independent of clustering) and write again
 	Point<float> bl(0.0, 0.0);
 	Point<float> tr(1.0, 1.0);

 	FruchtermanReingold fdLayouter(bl, tr);
 	fdLayouter.draw(G);
 	PostscriptWriter psWriter2(G, true);
 	psWriter2.write("output/testForceGraph.eps");
}

 TEST_F(VizGTest, tryMaxentLayouter) {
  	// create graph
  	count n = 120;
  	count numClusters = 3;
  	double pin = 0.175;
  	double pout = 0.005;

  	GraphGenerator graphGen;
  	Graph G = graphGen.makeClusteredRandomGraph(n, numClusters, pin, pout);
  	G.initCoordinates();
  	INFO("Number of edges: ", G.numberOfEdges());

  	// draw (independent of clustering) and write again
  	Point<float> bl(0.0, 0.0);
  	Point<float> tr(1.0, 1.0);

  	MaxentStress msLayouter(bl, tr);
 	msLayouter.draw(G);
  	PostscriptWriter psWriter3(G, true);
  	psWriter3.write("output/testMaxentGraph.eps");
}

 TEST_F(VizGTest, testMultilevelLayouter) {
  	// create graph
  	count n = 120;
  	count numClusters = 3;
  	double pin = 0.175;
  	double pout = 0.005;

  	GraphGenerator graphGen;
  	Graph G = graphGen.makeClusteredRandomGraph(n, numClusters, pin, pout);
  	G.initCoordinates();
  	INFO("Number of edges: ", G.numberOfEdges());

  	// draw (independent of clustering) and write again
  	Point<float> bl(0.0, 0.0);
  	Point<float> tr(1.0, 1.0);

  	MultilevelLayouter mlLayouter(bl, tr);
  	mlLayouter.draw(G);
  	PostscriptWriter psWriter4(G, true);
  	psWriter4.write("output/testMultilevelGraph.eps");
 }

 TEST_F(VizGTest, testAllLayouters) {
 	// create graph
 	count n = 120;
 	count numClusters = 3;
 	double pin = 0.175;
 	double pout = 0.005;

 	GraphGenerator graphGen;
 	Graph G = graphGen.makeClusteredRandomGraph(n, numClusters, pin, pout);
 	G.initCoordinates();
 	INFO("Number of edges: ", G.numberOfEdges());

 	// draw (independent of clustering) and write again
 	Point<float> bl(0.0, 0.0);
 	Point<float> tr(1.0, 1.0);

 	FruchtermanReingold fdLayouter(bl, tr);
 	fdLayouter.draw(G);
 	PostscriptWriter psWriter2(G, true);
 	psWriter2.write("output/testForceGraph.eps");

  	MultilevelLayouter mlLayouter(bl, tr);
  	mlLayouter.draw(G);
  	PostscriptWriter psWriter4(G, true);
  	psWriter4.write("output/testMultilevelGraph.eps");
}

 TEST_F(VizGTest, testMultilevelDrawing) {
  	// read graph
	METISGraphReader reader;
	Graph G = reader.read("input/jazz.graph");

  	// draw
  	Point<float> bl(0.0, 0.0);
  	Point<float> tr(1.0, 1.0);
   	MultilevelLayouter mlLayouter(bl, tr);
   	mlLayouter.draw(G);
   	PostscriptWriter psWriter4(G, true);
   	psWriter4.write("output/testJazzMl.eps");
 }




} /* namespace NetworKit */


#endif /*NOGTEST */


