/*
 * MapperGTest.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#include "MapperGTest.h"

#include "../StaticMapper.h"
#include "../RecBisMapper.h"
#include "../RcmMapper.h"
#include "../RcmMapperWW.h"
#include "../GreedyMapper.h"

#include "../../community/GraphClusteringTools.h"
#include "../../partitioning/BalancedLabelPropagation.h"
#include "../../graph/Graph.h"
#include "../../graph/GraphGenerator.h"
#include "../../io/DibapGraphReader.h"
#include "../../io/METISGraphReader.h"
#include "../../io/METISGraphWriter.h"
#include "../../community/ClusteringGenerator.h"






namespace NetworKit {


MapperGTest::MapperGTest() {

}

MapperGTest::~MapperGTest() {

}


#if 0
TEST_F(MapperGTest, tryRcmMapping) {
	// read application graph
	METISGraphReader graphReader;
	Graph appGraph = graphReader.read("input/wing.graph");
	count k = 25;

	// generate or read community/partition
	BalancedLabelPropagation partitioner(1.75);
	Partition partition = partitioner.run(appGraph, k);
//	ClusteringGenerator clusteringGenerator;
//	Partition partition = clusteringGenerator.makeContinuousBalancedClustering(appGraph, k);

	// read host (processor) graph
	Graph host = graphReader.read("input/mapping/grid-5x5-dist-arch.graph");

	// compute communication graph
	Graph commGraph = GraphClusteringTools::communicationGraph(appGraph, partition);

	// evaluate trivial mapping
	RcmMapper mapper;
	Mapping mapping = mapper.trivial(commGraph, host);
	edgeweight cost = mapper.cost(commGraph, host, mapping);
	INFO("Cost of trivial mapping (airfoil1 " , k , " parts onto 5x5 grid): " , cost);

	// call RCM mapping routine
	mapping = mapper.run(commGraph, host);

	// check and evaluate mapping
	cost = mapper.cost(commGraph, host, mapping);
	INFO("Cost of RCM mapping (airfoil1 " , k , " parts onto 5x5 grid): " , cost);
}
#endif 

TEST_F(MapperGTest, tryCommunicationGraph) {
	// TODO: read graph
	Graph g;

	// TODO: generate or read community/partition
	Partition partition;

	// TODO: compute communication graph
	Graph commGraph = GraphClusteringTools::communicationGraph(g, partition);

	// TODO: check communication graph

}


} // namespace EnsembleClustering

#endif
