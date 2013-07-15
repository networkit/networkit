/*
 * MapperGTest.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#include "MapperGTest.h"


namespace NetworKit {


MapperGTest::MapperGTest() {

}

MapperGTest::~MapperGTest() {

}


TEST_F(MapperGTest, tryRcmMapping) {
	// TODO: read graph
	Graph g;

	// TODO: read host (processor) graph
	Graph host;

	// TODO: generate or read clustering/partition
	Clustering partition;

	// TODO: compute communication graph
	Graph commGraph = partition.communicationGraph(g);

	// call mapping routine
	RcmMapper mapper;
	std::map<index, index> mapping = mapper.run(commGraph, host);

	// TODO: check and evaluate mapping
}


TEST_F(MapperGTest, tryCommunicationGraph) {
	// TODO: read graph
	Graph g;

	// TODO: generate or read clustering/partition
	Clustering partition;

	// TODO: compute communication graph
	Graph commGraph = partition.communicationGraph(g);

	// TODO: check communication graph

}


} // namespace EnsembleClustering

#endif
