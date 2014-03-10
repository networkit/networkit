/*
 * GraphProperties.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "GraphProperties.h"


namespace NetworKit {

GraphProperties::GraphProperties() {

}

GraphProperties::~GraphProperties() {

}

std::vector<count> GraphProperties::degreeDistribution(Graph& G) {
	count maxDegree = minMaxDegree(G).second;
	std::vector<count> distribution(maxDegree+1, 0);
	G.forNodes([&](node v){
		count i = G.degree(v);
		distribution[i]++;
	});
	return distribution;
}


std::vector<double> GraphProperties::localClusteringCoefficients(Graph& G) {
	count n = G.numberOfNodes();
	std::vector<double> numerator(n); //
	std::vector<double> denominator(n); // $\deg(u) \cdot ( \deg(u) - 1 )$
	std::vector<double> coefficient(n); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

	G.balancedParallelForNodes([&](node u){
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

	G.balancedParallelForNodes([&](node u){
		denominator[u] = G.degree(u) * (G.degree(u) - 1);
	});

	G.balancedParallelForNodes([&](node u){
		if (denominator[u] == 0.0) {
			coefficient[u] = 0.0;
		} else {
			coefficient[u] = numerator[u] / denominator[u];
		}
	});

	return coefficient;
}

std::vector<double> GraphProperties::localClusteringCoefficientPerDegree(Graph& G) {

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

double GraphProperties::averageLocalClusteringCoefficient(Graph& G) {
	std::vector<double> coefficients = GraphProperties::localClusteringCoefficients(G);
	double sum = 0.0;
	for (double c : coefficients) {
		sum += c;
	}
	double avg = sum / G.numberOfNodes();
	return avg;
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

std::vector<unsigned int> GraphProperties::degreeSequence(Graph& G) {
	std::vector<unsigned int> sequence(G.numberOfNodes()); // TODO: revert to count when cython issue fixed

	G.parallelForNodes([&](node v) {
		sequence[v] = G.degree(v);
	});

	return sequence;
}

double GraphProperties::averageDegree(const Graph& G) {

	count n = G.numberOfNodes();

	count degSum = G.parallelSumForNodes([&](node v){
		return G.degree(v);
	});

	double avgDeg = degSum / (double) n;
	return avgDeg;
}

double GraphProperties::degreeAssortativity(const Graph& G, bool useWeights) {
	// note: a parallel implementation would rather follow Newman's book, p. 267

	double r = 0.0; // result
	double A = 0.0; // accumulates degree products
	double B = 0.0; // accumulates degree sums
	double C = 0.0; // accumulates sum of degree squares

	double degu = 0.0; // temp storage for degree(u)
	double degv = 0.0; // temp storage for degree(v)

	double halfVolume = 0.0; // if needed, halfVolume accumulates the total edge weight of the graph (such a routine exists, but is not called for performance reasons)

	// iterate over edges and accumulate
	if (G.isMarkedAsWeighted() && useWeights) {
		G.forWeightedEdges([&](node u, node v, edgeweight ew) {
			degu = G.weightedDegree(u);
			degv = G.weightedDegree(v);
			A += degu * degv;
			B += degu + degv;
			C += degu*degu + degv*degv;
			halfVolume += ew;
		});
	}
	else {
		G.forEdges([&](node u, node v) {
			degu = G.degree(u);
			degv = G.degree(v);
			A += degu * degv;
			B += degu + degv;
			C += degu*degu + degv*degv;
		});

		halfVolume = G.numberOfEdges();
	}

	double volume = 2.0 * halfVolume;
	A = A / halfVolume;
	B = B / volume;
	B = B*B;
	C = C / volume;

	TRACE("A: ", A, ", B: ", B, ", C: ", C);

	assert(C != B);
	r = (A - B) / (C - B);
	return r;
}


} /* namespace NetworKit */
