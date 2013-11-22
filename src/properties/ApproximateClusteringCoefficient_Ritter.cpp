/*
 * ApproximateClusteringCoefficient_Ritter.cpp
 */

#include <stdlib.h>
#include <time.h>
#include <random>

#include "ApproximateClusteringCoefficient_Ritter.h"

namespace NetworKit {

ApproximateClusteringCoefficient_Ritter::ApproximateClusteringCoefficient_Ritter() {
	srand(time(NULL));
}

ApproximateClusteringCoefficient_Ritter::~ApproximateClusteringCoefficient_Ritter() {
}

double ApproximateClusteringCoefficient_Ritter::calculate(Graph& G, count k) {
	const count n = G.numberOfNodes();
	count degree[n];
	count cumulated_propabilities[n];
	count propability_sum = 0L;

	G.forNodes([&](node v) {
		degree[v] = G.degree(v);
		if (degree[v] > 1) {
			propability_sum += degree[v] * (degree[v] - 1) / 2;
		}
		cumulated_propabilities[v] = propability_sum;
	});

	count triangles_found = 0;

	for (index i = 1; i <= k; i++) {
		node v = randomNode(cumulated_propabilities, n, propability_sum);
		node w = G.randomNeighbor(v);
		node u;
		do {
			u = G.randomNeighbor(v);
		} while (u == w);

		if (G.hasEdge(w, u)) {
			triangles_found++;
		}

		if (i % 100000 == 0) {
			std::cout << "i = " << i << ", coefficient = " << (double) triangles_found / i << "\n";
		}
	}

	return (double) triangles_found / k;
}

/*
 * Binary search for the index i that fulfils array[i - 1] < key <= array[i]
 */
index ApproximateClusteringCoefficient_Ritter::indexOfBiggestValueLoEq(count* array, count size, count key) {
	// edge cases
	if (key <= array[0])
		return 0;
	if (key > array[size - 2])
		return size - 1;

	index left = 0;
	index right = size - 1;

	while (left <= right) {
		index middle = (left + right) / 2;
		if (key > array[middle]) {
			left = middle + 1;
		} else {
			// key <= array[middle]
			if (middle == 0 || key > array[middle-1]) {
				// array[middle - 1] > key => array[middle]
				return middle;
			}
			right = middle - 1;
		}
	}

	return -1;
}

count uniformRandom(count min, count max) {
	static count offset = 0;

	count currentMax = 1;
	count currentValue = 0;
	while(currentMax < max) {
		currentValue = currentValue * RAND_MAX + rand();
		currentMax *= RAND_MAX;
	}

	count value = currentValue % max;
	return offset = (value + offset) % max;
}

index ApproximateClusteringCoefficient_Ritter::randomNode(count* cumulated_propabilities, count n, count propability_sum) {
	// std::random_device rd;
	// std::mt19937 gen(rd());
	// return std::uniform_int_distribution<count>(0, propability_sum);
	return indexOfBiggestValueLoEq(cumulated_propabilities, n, uniformRandom(0, propability_sum));
}


} /* namespace NetworKit */