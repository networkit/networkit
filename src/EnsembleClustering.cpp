//============================================================================
// Name        : EnsembleClustering.cpp
// Author      : Christian Staudt
// Version     :
// Copyright   : © 2012, Christian Staudt
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <utility>
#include <unordered_map>

// log4cxx
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"


// cppunit
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>

// GoogleTest
#include "gtest/gtest.h"

// EnsembleClustering
#include "aux/log.h"
#include "aux/Noise.h"
#include "graph/Graph.h"
#include "input/METISParser.h"
#include "input/METIStoSTINGER.h"
#include "matching/Matching.h"
#include "clustering/Clustering.h"
#include "clustering/ClusteringGenerator.h"
#include "graph/GraphGenerator.h"
#include "clustering/Modularity.h"

extern "C" {
#include "stinger.h"
}




using namespace EnsembleClustering;




void testMETIStoSTINGER() {

	std::string graphPath = "/Users/cls/workspace/Data/DIMACS/kron_g500-simple-logn16.graph";
	LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), "trying to read from graph file " << graphPath);


	Graph* G;
	METIStoSTINGER* m2s = new METIStoSTINGER();
	G = m2s->read(graphPath);

	LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), "read graph " << G << " from file " << graphPath);

}




void testMatching() {
	INFO("testing matching");

	int n = 10e7;
	Matching M(n);

	#pragma omp parallel for
	for (node u = 1; u <= n; ++u) {
		M.match(u, (u + 1) % n);
	}

//	std::cout << "Node " << u << " is matched: " << M.isMatched(u) << std::endl;
//	std::cout << "Node " << u << " is matched with " << M[u] << std::endl;

}




/**
 * Make a complete graph with n vertices.
 *
 */
Graph& makeCompleteGraph(int n) {

	Graph G;

	for (node u = 0; u < n; ++u) {
		for (node v = u + 1; v < n; ++v) {
			G.insertEdge(u, v);
		}
	}

	DEBUG("number of edges " << G.numberOfEdges());

	return G;
}


/**
 * Call this first to configure logging output.
 */
void configureLogging() {
	// configure logging
	log4cxx::BasicConfigurator::configure();
	log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getDebug());
}


/**
 * Call this to run all cppunit unit tests.
 */
int runUnitTestsCppUnit() {
	// Informiert Test-Listener ueber Testresultate
	CPPUNIT_NS::TestResult testresult;

	// Listener zum Sammeln der Testergebnisse registrieren
	CPPUNIT_NS::TestResultCollector collectedresults;
	testresult.addListener(&collectedresults);

	// Listener zur Ausgabe der Ergebnisse einzelner Tests
	CPPUNIT_NS::BriefTestProgressListener progress;
	testresult.addListener(&progress);

	// Test-Suite ueber die Registry im Test-Runner einfuegen
	CPPUNIT_NS::TestRunner testrunner;
	testrunner.addTest(CPPUNIT_NS::TestFactoryRegistry::getRegistry().makeTest());
	testrunner.run(testresult);

	// Resultate im Compiler-Format ausgeben
	CPPUNIT_NS::CompilerOutputter compileroutputter(&collectedresults, std::cerr);
	compileroutputter.write();

	// Rueckmeldung, ob Tests erfolgreich waren
	return collectedresults.wasSuccessful() ? 0 : 1;
}

/**
 * Call this to run all GoogleTest unit tests
 */
int runUnitTests(int argc, char **argv) {
::testing::InitGoogleTest(&argc, argv);
 return RUN_ALL_TESTS();
}





void testModularity() {
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
}





int main(int argc, char **argv) {

	std::cout << "running EnsembleClustering" << std::endl;

	configureLogging();

	 ::testing::InitGoogleTest(&argc, argv);
	 return RUN_ALL_TESTS();


	return 0;
}
