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


std::pair<count, count> GraphProperties::estimatedDiameterRange_Hoske(const Graph& G, double error) {
	using namespace std;

	/* Determines a node v with dist(u, v) = ecc_G(u) and the dist(u, v). */
	auto compute_ecc = [&] (const Graph& G, node u) -> pair<node, count> {
		static BFS bfs;
		auto dists = bfs.run(G, u);
		auto max_iter = max_element(begin(dists), end(dists));
		return {distance(begin(dists), max_iter), *max_iter};
	};

	/* BFS that calls f with the visited edges. */
	/* Note: the function Graph::breadthFirstEdgesFrom that should
	   do the same has not been implemented! */
	auto bfs_edges = [&] (const Graph& G, node u, function<void(node, node)> f) {
		queue<node> q;
		vector<bool> visited(G.numberOfNodes(), false);
		q.push(u);
		visited[u] = true;

		while (!q.empty()) {
			node x = q.front(); q.pop();
			G.forNeighborsOf(x, [&] (node y) {
				if (!visited[y]) {
					f(x, y);
					visited[y] = true;
					q.push(y);
				}
			});
		}
	};

	/* Diameter estimate: lowerBounds <= diam(G) <= upperBound. */
	count lowerBound = 0;
	count upperBound = numeric_limits<count>::max();
	const count n = G.numberOfNodes();

	/* Nodes sorted by decreasingly by degree. */
	vector<node> high_deg(n);
	iota(begin(high_deg), end(high_deg), 0);
	sort(begin(high_deg), end(high_deg), [&] (node u, node v) {
		return G.degree(u) > G.degree(v);
	});

	/* Random node. */
	static const default_random_engine random;
	auto random_node = bind(uniform_int_distribution<node>(0, n - 1), random);

	/* While not converged: update estimate. */
	count niter = 0;
	while ((upperBound - lowerBound) >= error*lowerBound && niter < n) {
		count ecc;

		/* ecc(u) <= diam(G) */
		node u = random_node();
		tie(ignore, ecc) = compute_ecc(G, u);
		lowerBound = max(lowerBound, ecc);

		/* diam(G) <= diam(BFS_Tree(v)) */
		node v = high_deg[niter];
		Graph bfs_tree(n);
		bfs_edges(G, v, [&] (node a, node b) {
			bfs_tree.addEdge(a, b);
		});
		/* Note: A bit inefficient. One BFS per iteration too many. */
		node w;
		tie(w, ecc) = compute_ecc(G, u);
		lowerBound = max(lowerBound, ecc);
		/* diam(T) = ecc_T(w) by problem 4. */
		tie(ignore, ecc) = compute_ecc(bfs_tree, w);
		upperBound = min(upperBound, ecc);

		niter++;
	}

	return {lowerBound, upperBound};
}

} /* namespace NetworKit */

