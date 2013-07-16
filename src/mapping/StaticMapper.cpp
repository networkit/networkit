/*
 * StaticMapper.cpp
 *
 *  Created on: Jul 2, 2013
 *      Author: Henning
 */

#include "StaticMapper.h"

namespace NetworKit {

StaticMapper::StaticMapper() {

}

StaticMapper::~StaticMapper() {

}

Mapping StaticMapper::trivial(Graph& guest, Graph& host) {
	Mapping mapping;
	assert(guest.numberOfNodes() <= host.numberOfNodes());

	guest.forNodes([&](node v) {
		mapping.insert(std::make_pair(v, v));
	});

	return mapping;
}

edgeweight StaticMapper::cost(const Graph& guest, const Graph& host, Mapping& mapping) {
	edgeweight cost = 0.0;

	guest.forWeightedEdges([&](node u, node v, edgeweight w) {
		cost += w * host.weight(mapping[u], mapping[v]); // TODO: assumes complete graph, otherwise use distance!
	});

	return cost;
}

} /* namespace NetworKit */
