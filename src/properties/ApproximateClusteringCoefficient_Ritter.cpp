/*
 * ApproximateClusteringCoefficient_Ritter.cpp
 */

#include <stdlib.h>
#include <random>

#include "ApproximateClusteringCoefficient_Ritter.h"

namespace NetworKit {

ApproximateClusteringCoefficient_Ritter::ApproximateClusteringCoefficient_Ritter() {
}

ApproximateClusteringCoefficient_Ritter::~ApproximateClusteringCoefficient_Ritter() {
}

double ApproximateClusteringCoefficient_Ritter::calculate(const Graph& G, count k) {
	const count n = G.numberOfNodes();
	std::vector<count> degree(n);
	std::vector<count> cumulatedPropabilities(n);
	count propabilitySum = 0L;

	G.forNodes([&](node v) {
		degree[v] = G.degree(v);
		if (degree[v] > 1) {
			propabilitySum += degree[v] * (degree[v] - 1) / 2;
		}
		cumulatedPropabilities[v] = propabilitySum;
		cumulatedPropabilities[v] = degree[v] * (degree[v] - 1) / 2;
	});

	std::mt19937_64 randGen;
	// std::uniform_int_distribution<node> propabilityDistribution(0, propabilitySum - 1);
	std::discrete_distribution<count> propabilityDistribution(cumulatedPropabilities.begin(), cumulatedPropabilities.end());

	count triangles_found = 0;

	for (index i = 1; i <= k; i++) {
		// node v = indexOfBiggestValueLoEq(cumulatedPropabilities, propabilityDistribution(randGen));
		node v = propabilityDistribution(randGen);
		node w = G.randomNeighbor(v);
		node u;
		do {
			u = G.randomNeighbor(v);
		} while (u == w);

		if (G.hasEdge(w, u)) {
			triangles_found++;
		}

		// if (i % 10000 == 0) {
		// 	std::cout << "i = " << i << ", coefficient = " << (double) triangles_found / i << "\n";
		// }
	}

	return (double) triangles_found / k;
}

/*
 * Binary search for the index i that fulfills array[i - 1] < key <= array[i]
 */
index ApproximateClusteringCoefficient_Ritter::indexOfBiggestValueLoEq(std::vector<count> values, count key) {
	// edge cases
	if (key <= values[0])
		return 0;
	if (key > values[values.size() - 2])
		return values.size() - 1;

	index left = 0;
	index right = values.size() - 1;

	while (left <= right) {
		index middle = (left + right) / 2;
		if (key > values[middle]) {
			left = middle + 1;
		} else {
			// key <= values[middle]
			if (middle == 0 || key > values[middle-1]) {
				// values[middle - 1] > key => values[middle]
				return middle;
			}
			right = middle - 1;
		}
	}

	return -1;
}

} /* namespace NetworKit */