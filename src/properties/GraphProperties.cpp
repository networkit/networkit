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
	std::vector<count> distribution(maxDegree, 0);
	G.forNodes([&](node v){
		count i = G.degree(v);
		distribution[i]++;
	});
}


std::vector<double> GraphProperties::localClusteringCoefficients(Graph& G) {
	count n = G.numberOfNodes();
	std::vector<count> numerator(n); //
	std::vector<count> denominator(n); // $\deg(u) \cdot ( \deg(u) - 1 )$
	std::vector<double> coefficient(n); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$


	G.parallelForNodes([&](node u){
		count edgeCount = 0;
		G.forEdgesOf(u, [&](node u, node v) {
			G.forEdgesOf(v, [&](node v, node w){
				if (G.hasEdge(u, w)) {
					edgeCount += 1;
				}
			});
		});

		numerator[u] = edgeCount; // factor 2 is omitted because each edge has been counted twice
	});

	G.parallelForNodes([&](node u){
		denominator[u] = G.degree(u) * (G.degree(u) - 1);
	});

	G.parallelForNodes([&](node u){
		coefficient[u] = numerator[u] / denominator[u];
	});

}

std::vector<double> GraphProperties::localClusteringCoefficientPerDegree(Graph& G) {

	std::vector<count> degDist = degreeDistribution(G);
	std::vector<double> coefficient = localClusteringCoefficients(G);

	std::vector<double> perDegree(degDist.size(), 0.0);

	G.forNodes([&](node u){
		perDegree[G.degree(u)] += coefficient[u];
	});

	// get the average local clustering coefficient for nodes of each degreee
	for (index i = 0; i < degDist.size(); ++i) {
		perDegree[i] = perDegree[i] / (double) degDist[i];
	}

	return perDegree;
}



std::pair<count, count> GraphProperties::minMaxDegree(Graph& G) {

	count min = G.numberOfNodes();
	count max = 0;

	G.forNodes([&](node v){
		count d = G.degree(v);
		if (d < min) {
			min = d;
		}
		if (d > max) {
			max = d;
		}
	});

	return std::pair<count, count>(min, max);
}

} /* namespace NetworKit */
