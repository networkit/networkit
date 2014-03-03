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
	Partition zeta = clusteringGen.makeRandomClustering(G, numClusters);

	// create coordinates
	G.forNodes([&](node u) {
		Point<float> p(drand48(), drand48());
		G.setCoordinate(u, p);
	});

	// write graph to file
	PostscriptWriter psWriter(G);
	psWriter.write(zeta, "testGraph.eps");
}

static float edgeDistanceSum(Graph& G) {
	float dist = 0.0f;

	G.forEdges([&](node u, node v) {
		Point<float> p = G.getCoordinate(u) - G.getCoordinate(v);
		dist += p.length();
	});

	return dist;
}

TEST_F(VizGTest, testFRLayouter) {
 	// create graph
 	count n = 80;
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

 	// test edge distances
 	float dist = edgeDistanceSum(G);
 	float avg = dist / (float) G.numberOfEdges();
 	INFO("avg edge length: ", avg);
 	EXPECT_LE(avg, 0.25);
}

 TEST_F(VizGTest, tryMaxentLayouter) {
  	// create graph
  	count n = 80;
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

 	// test edge distances
 	float dist = edgeDistanceSum(G);
 	float avg = dist / (float) G.numberOfEdges();
 	DEBUG("avg edge length: ", avg);
 	EXPECT_LE(avg, 0.25);
}

 TEST_F(VizGTest, testMultilevelLayouter) {
  	// create graph
  	count n = 80;
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

  	// test edge distances
 	float dist = edgeDistanceSum(G);
 	float avg = dist / (float) G.numberOfEdges();
 	DEBUG("avg edge length: ", avg);
 	EXPECT_LE(avg, 0.25);
}

 TEST_F(VizGTest, testGraphDrawing) {
 	// create graph
	METISGraphReader reader;
	Graph G = reader.read("input/lesmis.graph");

 	// draw (independent of clustering) and write again
 	Point<float> bl(0.0, 0.0);
 	Point<float> tr(1.0, 1.0);

 	FruchtermanReingold fdLayouter(bl, tr);
 	fdLayouter.draw(G);
 	PostscriptWriter psWriter2(G, true);
 	psWriter2.write("output/testLesmisFR.eps");

 	// test edge distances
 	float dist = edgeDistanceSum(G);
 	float avg = dist / (float) G.numberOfEdges();
 	INFO("avg edge length: ", avg);
 	EXPECT_LE(avg, 0.25);

 	MultilevelLayouter mlLayouter(bl, tr);
  	mlLayouter.draw(G);
  	PostscriptWriter psWriter4(G, true);
  	psWriter4.write("output/testLesmisMl.eps");

 	// test edge distances
 	dist = edgeDistanceSum(G);
 	avg = dist / (float) G.numberOfEdges();
 	INFO("avg edge length: ", avg);
 	EXPECT_LE(avg, 0.25);
}


} /* namespace NetworKit */


#endif /*NOGTEST */


