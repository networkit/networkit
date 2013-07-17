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
	// read application graph
	METISGraphReader graphReader;
	Graph appGraph = graphReader.read("input/airfoil1.graph");
	count k = 512;

	// generate or read clustering/partition
	ClusteringGenerator clusteringGenerator;
	Clustering partition = clusteringGenerator.makeContinuousBalancedClustering(appGraph, k);

	// read host (processor) graph
	Graph host = graphReader.read("input/mapping/grid-8x8x8-dist-arch.graph");

	// compute communication graph
	Graph commGraph = partition.communicationGraph(appGraph);

	// evaluate trivial mapping
	RcmMapperWW mapper;
	Mapping mapping = mapper.trivial(commGraph, host);
	edgeweight cost = mapper.cost(commGraph, host, mapping);
	INFO("Cost of RCM mapping (airfoil1 512 parts onto 8x8x8 grid): " << cost);

	// call RCM mapping routine
	mapping = mapper.run(commGraph, host);

	// check and evaluate mapping
	cost = mapper.cost(commGraph, host, mapping);
	INFO("Cost of RCM mapping (airfoil1 512 parts onto 8x8x8 grid): " << cost);
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
