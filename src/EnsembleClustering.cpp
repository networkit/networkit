//============================================================================
// Name        : EnsembleClustering.cpp
// Author      : Christian Staudt
// Version     :
// Copyright   : © 2012, Christian Staudt
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <utility>

#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"

#include "aux/log.h"
#include "aux/Noise.h"
#include "graph/Graph.h"
#include "input/METISParser.h"
#include "input/METIStoSTINGER.h"
#include "matching/Matching.h"

extern "C" {
#include "stinger.h"
}


using namespace EnsembleClustering;



void testMETISParser() {

	METISParser* parser = new METISParser();

	parser->open("/Users/cls/workspace/Data/DIMACS/kron_g500-simple-logn16.graph");
	std::pair<int, int> header = parser->getHeader();

	int lc = 0;
	while (parser->hasNext()) {
		parser->getNext();
		lc++;
	}

	LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), "parsed " << lc << " lines");

	parser->close();
}


void testMETIStoSTINGER() {

	std::string graphPath = "/Users/cls/workspace/Data/DIMACS/kron_g500-simple-logn16.graph";
	LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), "trying to read from graph file " << graphPath);


	Graph* G;
	METIStoSTINGER* m2s = new METIStoSTINGER();
	G = m2s->read(graphPath);

	LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), "read graph " << G << " from file " << graphPath);

}


void testNoise() {

	INFO("testing noise");
	Noise noise(-0.5, 0.5);
	double x = 1.0;
	for (int i = 0; i < 10; i++) {
		std::cout << noise.add(x) << std::endl;
	}
}


void testMatching() {
	INFO("testing matching");

	int n = 10e7;
	Matching M(n);

	#pragma omp parallel for
	for (node u = 0; u < n; ++u) {
		M.match(u, (u + 1) % n);
	}

//	std::cout << "Node " << u << " is matched: " << M.isMatched(u) << std::endl;
//	std::cout << "Node " << u << " is matched with " << M[u] << std::endl;

}


/**
 * Test iteration
 */
void testIteration(Graph& G) {
	INFO("testing iteration");

	int64_t etype = G.defaultEdgeType;

	STINGER_PARALLEL_FORALL_EDGES_BEGIN(G.asSTINGER(), etype) {
		node u = STINGER_EDGE_SOURCE;
		node v = STINGER_EDGE_DEST;
		std::printf("found edge (%d, %d) with weight %f \n", u, v, stinger_edgeweight(G.asSTINGER(), u, v, etype));
	} STINGER_PARALLEL_FORALL_EDGES_END();

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


int main() {

	std::cout << "running EnsembleClustering" << std::endl;

	// configure logging
	log4cxx::BasicConfigurator::configure();
	log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getDebug());

	INFO("test debug macro");

	int n = 10;
	Graph G = makeCompleteGraph(n);
	testIteration(G);

	// testMatching();


	return 0;
}
