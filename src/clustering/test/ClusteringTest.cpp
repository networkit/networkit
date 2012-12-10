/*
 * ClusteringTest.cpp
 *
 *  Created on: 06.12.2012
 *      Author: cls
 */

#include "ClusteringTest.h"

namespace EnsembleClustering {

CPPUNIT_TEST_SUITE_REGISTRATION(ClusteringTest);

ClusteringTest::ClusteringTest() {
	// TODO Auto-generated constructor stub

}

ClusteringTest::~ClusteringTest() {
	// TODO Auto-generated destructor stub
}

void ClusteringTest::setUp() {
	this->randomGraph = this->gen.makeErdosRenyiGraph(20, 0.2);
}

void ClusteringTest::tearDown() {
}

void ClusteringTest::testModularity() {

	GraphGenerator graphGenerator;

	Graph G = graphGenerator.makeErdosRenyiGraph(100, 0.1);

	ClusteringGenerator clusteringGenerator;

	Clustering singleton = clusteringGenerator.makeSingletonClustering(G);
	Clustering one = clusteringGenerator.makeOneClustering(G);

	Modularity modularity(G);

	double modSingleton = modularity.getQuality(singleton);
	double modOne = modularity.getQuality(one);

	CPPUNIT_ASSERT_EQUAL_MESSAGE("A singleton clustering should have a modularity of 0", modSingleton, 0.0);
	CPPUNIT_ASSERT_EQUAL_MESSAGE("A 1-clustering should have a modularity of 0", modOne, 0.0);
}

} /* namespace EnsembleClustering */
