/*
 * GraphProperties.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "GraphProperties.h"
#include "GraphProperties_Ritter.h"

#include <list>
#include <random>

namespace NetworKit {

GraphProperties::GraphProperties() {

}

GraphProperties::~GraphProperties() {

}

std::vector<count> GraphProperties::degreeDistribution(const Graph& G) {
	count maxDegree = minMaxDegree(G).second;
	std::vector<count> distribution(maxDegree+1, 0);
	G.forNodes([&](node v){
		count i = G.degree(v);
		distribution[i]++;
	});
	return distribution;
}


std::vector<double> GraphProperties::localClusteringCoefficients(const Graph& G) {
	count n = G.numberOfNodes();
	std::vector<double> numerator(n); //
	std::vector<double> denominator(n); // $\deg(u) \cdot ( \deg(u) - 1 )$
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
		if (denominator[u] == 0.0) {
			coefficient[u] = 0.0;
		} else {
			coefficient[u] = numerator[u] / denominator[u];
		}
	});

	return coefficient;
}

std::vector<double> GraphProperties::localClusteringCoefficientPerDegree(const Graph& G) {

	std::vector<count> degDist = degreeDistribution(G);
	std::vector<double> coefficient;
	std::vector<double> perDegree(degDist.size(), 0.0);

	if (G.numberOfNodes() > 1 ) {
		coefficient = localClusteringCoefficients(G);

			G.forNodes([&](node u){
				perDegree[G.degree(u)] += coefficient[u];
			});

			// get the average local clustering coefficient for nodes of each degreee
			for (index i = 2; i < degDist.size(); ++i) {
				if (degDist[i] == 0) {
					perDegree[i] = 0.0; // TODO: should this be -1
				} else {
					perDegree[i] = perDegree[i] / (double) degDist[i];
				}
			}
	}



	// allows to avoid the situation, when local clustering coefficient is calculated for 0-1 degree nodes.
	// These nodes are warranted not to be triangle centers, thus we avoid calculating the coefficients for the,
	degDist[0] = 0;
	if (G.numberOfNodes() > 0 ) degDist[1] = 0;

	return perDegree;
}

double GraphProperties::averageLocalClusteringCoefficient(const Graph& G) {
	std::vector<double> coefficients = GraphProperties::localClusteringCoefficients(G);
	double sum = 0.0;
	for (double c : coefficients) {
		sum += c;
	}
	double avg = sum / G.numberOfNodes();
	return avg;
}

std::pair<count, count> GraphProperties::minMaxDegree(const Graph& G) {

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


double GraphProperties::averageDegree(const Graph& G) {

	count n = G.numberOfNodes();

	count degSum = G.parallelSumForNodes([&](node v){
		return G.degree(v);
	});

	double avgDeg = degSum / (double) n;
	return avgDeg;
}

double GraphProperties::approximateGlobalClusteringCoefficient(
		const Graph& G, double approximationError, double probabilityError) {
	ApproximateClusteringCoefficient_Hoske acc;
	count numIters = acc.niters(approximationError, probabilityError);
	return acc.calculate(true, G, numIters);
}

std::pair<count, count> GraphProperties::estimatedDiameterRange_Ritter(const Graph& G) {
	static const count p = 5;

	const count n = G.numberOfNodes();
	std::mt19937_64 randGen;
	std::uniform_int_distribution<node> distribution(0, n - 1);

	if (G.numberOfEdges() == 0) {
		// empty graph
		return std::make_pair(0, 0);
	}

	if (GraphProperties_Ritter::eccentricity(G, 0) == GraphProperties_Ritter::INF_DIST) {
		// Graph G is not connected
		return std::make_pair(GraphProperties_Ritter::INF_DIST, GraphProperties_Ritter::INF_DIST);
	}

	count lowerBound = 1;
	count upperBound = std::numeric_limits<count>::max();

	const count maxDegree = minMaxDegree(G).second;

	std::list<node> nodesByDegree[maxDegree + 1];
	G.forNodes([&] (node v) {
		nodesByDegree[G.degree(v)].push_front(v);
	});

	index min_d = 0;
	index max_d = maxDegree;
	while (upperBound - lowerBound > p) {
		// let min_d and max_d point to non empty list again
		while (nodesByDegree[min_d].empty() && min_d < maxDegree) min_d++;
		while (nodesByDegree[max_d].empty() && max_d > 0) max_d--;

		// we need at least 2 nodes in nodesByDegree to continue
		if (max_d < min_d || (min_d == max_d && nodesByDegree[min_d].size() < 2)) {
			break;
		}

		// try to improve lower bound
		node v = distribution(randGen);

		// using the eccentricity of a low degree node v
		// node v = nodesByDegree[min_d].front();
		// nodesByDegree[min_d].pop_front();

		count ecc = GraphProperties_Ritter::eccentricity(G, v);
		if (ecc > lowerBound) {
			lowerBound = ecc;
		}

		// try to improve upper bound
		// using the diameter of a spanning tree T from a high degree node u
		node u = nodesByDegree[max_d].front();

		nodesByDegree[max_d].pop_front();
		Graph T = GraphProperties_Ritter::spanningTree(G, u);
		count dia = GraphProperties_Ritter::diameterOfTree(T, u);
		if (dia < upperBound) {
			upperBound = dia;
		}
		// printf("lower bound: %u, upper bound: %u\n", lowerBound, upperBound);
	}

	return std::make_pair(lowerBound, upperBound);
}

} /* namespace NetworKit */

