/*
 * Assortativity.cpp
 *
 *  Created on: Jun 13, 2015
 *      Author: Christian Staudt
 */

#include "Assortativity.h"

namespace NetworKit {

Assortativity::Assortativity(const Graph& G, const std::vector<double>& attribute) : G(G), emptyVector(), emptyPartition(), attribute(attribute),  partition(emptyPartition), nominal(false) {
	if (attribute.size() < G.upperNodeIdBound()) {
		throw std::runtime_error("attribute list has incorrect length: there must be an entry for each node");
	}
}


Assortativity::Assortativity(const Graph& G, const Partition& partition) : G(G), emptyVector(), emptyPartition(), partition(partition), attribute(emptyVector), nominal(true) {
	if (partition.numberOfElements() < G.upperNodeIdBound()) {
		throw std::runtime_error("partition has incorrect length: there must be an entry for each node");
	}
}


void Assortativity::run() {
	if (nominal) {
		// compact partition so matrix doesn't get unnecessarily large
		Partition P = partition;
		P.compact();
		// create kxk matrix with entries $e_{ij}$, the fraction of edges connecting nodes of type i to nodes of type j
		count k = P.upperBound();
		std::vector<std::vector<double>> E(k, std::vector<double>(k, 0.0));
		G.forEdges([&](node u, node v) {
			E[P[u]][P[v]] += 1;
		});
		// row and column sums $a_i$ and $b_i$
		std::vector<double> a(k, 0.0);
		std::vector<double> b(k, 0.0);
		// normalize and calculate sums
		count m = G.numberOfEdges();
		for (index i = 0; i < k; ++i) {
			for (index j = 0; j < k; ++j) {
				E[i][j] = E[i][j] / m;
				a[i] += E[i][j];
				b[j] += E[i][j];
			}
		}
		// calculate coefficient $r = \frac{\sum_i e_{ii} - \sum_i a_i b_i}{1 - \sum_i a_i b_i}$
		double diagSum = 0.0;
		double abSum = 0.0;
		for (index i = 0; i < k; ++i) {
			diagSum += E[i][i];
			abSum += a[i] * b[i];
		}
		double r = (diagSum - abSum) / (double) (1 - abSum);
		INFO("diagSum: ", diagSum);
		INFO("abSum: ", abSum);
		INFO("r: ", r);
		coefficient = r;
	} else {

	}

}


double Assortativity::getCoefficient() const {
	return coefficient;
}


}
