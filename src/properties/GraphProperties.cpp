/*
 * GraphProperties.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "GraphProperties.h"
#include "../graph/BFS.h"
#include "../graph/BFSTree.h"
#include <random>
#include <algorithm>


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


std::pair<count, count> GraphProperties::estimatedDiameterRange_Brueckner(
		const Graph& G) {

    using namespace std;

    count infDist = numeric_limits<count>::max();

    default_random_engine e;
    e.seed(random_device()());

    // Pick nodes randomly, weighted by their degree.
    vector<count> nodeWeights(G.numberOfNodes());
    G.forNodes([&](node v) {
            nodeWeights[v] = G.degree(v);
        });

	count lowerBound = 0;
	count upperBound = infDist;

    for (count i = 0;
         i < G.numberOfNodes() && upperBound - lowerBound > 5;
         i++) {

        discrete_distribution<node> nodeDistribution(nodeWeights.begin(), nodeWeights.end());
        node v = nodeDistribution(e);
        nodeWeights[v] = 0;
        BFSTree T(G, v);
        if (!T.spanning())
            return make_pair(infDist, infDist);
        else {
            lowerBound = max(lowerBound, T.depth());
            upperBound = min(upperBound, 2 * T.depth());
            BFS bfs;
            std::vector<count> dists = bfs.run(T, T.deepest());
            count maxDist = 0;
            for (auto dist : dists)
                maxDist = max(maxDist, dist);
            upperBound = min(upperBound, maxDist);
        }
    }

	return make_pair(lowerBound, upperBound);
}

count GraphProperties::exactDiameter_Brueckner(
        const Graph& G) {

    using namespace std;

    count diameter = 0;

	G.forNodesInRandomOrder([&](node v){
            BFS bfs;
            vector<count> distances = bfs.run(G, v);
            for (auto distance : distances) {
                if (diameter < distance)
                    diameter = distance;
            }
    });

    return diameter;
}

} /* namespace NetworKit */
