/*
 * BalancedPartitioner.cpp
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#include "BalancedPartitioner.h"

namespace NetworKit {

BalancedPartitioner::BalancedPartitioner(): balance(1.0) {

}

BalancedPartitioner::~BalancedPartitioner() {

}


Clustering BalancedPartitioner::multilevelRun(Graph& graph, count numParts) {
	DEBUG("Start level with " , graph.numberOfNodes() , " nodes");

	if (graph.numberOfNodes() <= 48 * numParts) {
		// terminate recursion
		return this->run(graph, numParts);
	}
	else {
		// contract
		MatchingContracter matchContract;
		LocalMaxMatcher matcher(none);
		Matching matching = matcher.run(graph);
		auto coarseInfo = matchContract.run(graph, matching);
		Graph& coarseGraph = coarseInfo.first;
		auto fineToCoarse = coarseInfo.second;

		// recurse
		Clustering coarsePartition = this->multilevelRun(coarseGraph, numParts);

		// prolongate, refine
		ClusteringProjector prolongator;
		Clustering partition = prolongator.projectBack(coarseGraph, graph,	fineToCoarse, coarsePartition);
		partition = this->rerun(graph, numParts, partition);

		// postsmooth quasi-deterministically and return
		return this->postsmooth(graph, numParts, partition);
	}
}

Clustering& BalancedPartitioner::multilevelRerun(Graph& graph, count numParts,
		Clustering& partition) {

	DEBUG("Start level with " , graph.numberOfNodes() , " nodes");

	if (graph.numberOfNodes() <= 48 * numParts) {
		// terminate recursion
		partition = this->rerun(graph, numParts, partition);
		return partition;
	}
	else {
		// contract
		MatchingContracter matchContract;
		LocalMaxMatcher matcher(none);
		Matching matching = matcher.run(graph);
		auto coarseInfo = matchContract.run(graph, matching);
		Graph& coarseGraph = coarseInfo.first;
		auto fineToCoarse = coarseInfo.second;

		// recurse
		Clustering coarsePartition(coarseGraph.numberOfNodes());
		graph.forNodes([&](node v) {
			coarsePartition[fineToCoarse[v]] = partition[v];
		});
		coarsePartition = this->multilevelRerun(coarseGraph, numParts, coarsePartition);

		// prolongate, refine
		ClusteringProjector prolongator;
		partition = prolongator.projectBack(coarseGraph, graph,	fineToCoarse, coarsePartition);
		partition = rerun(graph, numParts, partition);

		// postsmooth quasi-deterministically and return
		return postsmooth(graph, numParts, partition);
	}
}

void BalancedPartitioner::setBalance(float balanceFactor) {
	this->balance = balanceFactor;
}

} /* namespace NetworKit */
