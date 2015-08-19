/*
 * Assortativity.cpp
 *
 *  Created on: Jun 13, 2015
 *      Author: Christian Staudt
 */

#include "Assortativity.h"

namespace NetworKit {

Assortativity::Assortativity(const Graph& G, const std::vector<double>& attribute) : G(G), emptyVector(), emptyPartition(), attribute(attribute),  partition(emptyPartition) {
	if (attribute.size() < G.upperNodeIdBound()) {
		throw std::runtime_error("attribute list has incorrect length: there must be an entry for each node");
	}
}


Assortativity::Assortativity(const Graph& G, const Partition& partition) : G(G), emptyVector(), emptyPartition(), partition(partition), attribute(emptyVector) {
	if (partition.numberOfElements() < G.upperNodeIdBound()) {
		throw std::runtime_error("partition has incorrect length: there must be an entry for each node");
	}
}


void Assortativity::run() {
	throw std::runtime_error("TODO");
}


double Assortativity::getCoefficient() const {
	throw std::runtime_error("TODO");
	return 0.0;
}


}
