/*
 * ClusteringTest.cpp
 *
 *  Created on: 06.12.2012
 *      Author: cls
 */

#include "ClusteringTest.h"

namespace EnsembleClustering {

CPPUNIT_TEST_SUITE_REGISTRATION(ClusteringCppUnitTest);

ClusteringCppUnitTest::ClusteringCppUnitTest() {
	// TODO Auto-generated constructor stub

}

ClusteringCppUnitTest::~ClusteringCppUnitTest() {
	// TODO Auto-generated destructor stub
}

void ClusteringCppUnitTest::setUp() {
	this->randomGraph = this->gen.makeErdosRenyiGraph(20, 0.2);
}

void ClusteringCppUnitTest::tearDown() {
}

void ClusteringCppUnitTest::testModularity() {

	GraphGenerator graphGenerator;

	Graph G = graphGenerator.makeCompleteGraph(20);

	ClusteringGenerator clusteringGenerator;

	Clustering singleton = clusteringGenerator.makeSingletonClustering(G);
	Clustering one = clusteringGenerator.makeOneClustering(G);

	Modularity modularity(G);

	double modSingleton = modularity.getQuality(singleton);
	double modOne = modularity.getQuality(one);

	std::cout << "modSingleton: " << modSingleton << std::endl;
	std::cout << "modOne: " << modOne << std::endl;

	CPPUNIT_ASSERT_EQUAL(0.0, modSingleton);
	CPPUNIT_ASSERT_EQUAL(0.0, modOne);


//	CPPUNIT_ASSERT_EQUAL_MESSAGE("A singleton clustering should have a modularity of 0", 0.0, modSingleton);
//	CPPUNIT_ASSERT_EQUAL_MESSAGE("A 1-clustering should have a modularity of 0", 0.0, modOne);
}

} /* namespace EnsembleClustering */
