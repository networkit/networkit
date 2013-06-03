/*
 * GraphProperties.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "GraphProperties.h"

namespace NetworKit {

GraphProperties::GraphProperties() {
	// TODO Auto-generated constructor stub

}

GraphProperties::~GraphProperties() {
	// TODO Auto-generated destructor stub
}

std::vector<count> GraphProperties::degreeDistribution(Graph& G) {
	count maxDegree = minMaxDegree(G).second;
	vector<count> values (maxDegree, 0);

	G.forNodes([&](node v){
		count i = G.degree(v);
		vector[i]++;
	});

	return values;
}

std::vector<double> GraphProperties::localClusteringCoefficientPerDegree(
		Graph& G) {
}

} /* namespace NetworKit */
