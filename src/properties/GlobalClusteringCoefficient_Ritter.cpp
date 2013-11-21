/*
 * GlobalClusteringCoefficient_Ritter.cpp
 */

#include "GlobalClusteringCoefficient_Ritter.h"
#include <mutex>

namespace NetworKit {

GlobalClusteringCoefficient_Ritter::GlobalClusteringCoefficient_Ritter() {
}

GlobalClusteringCoefficient_Ritter::~GlobalClusteringCoefficient_Ritter() {
}

double GlobalClusteringCoefficient_Ritter::calculate(Graph& G) {
	std::mutex mtx;

	node numerator = 0;
	node denominator = 0;

	G.parallelForNodes([&](node v) {
		const node d = G.degree(v);
		if (d >= 2) {
			node n = 0;

			G.forTwoNeighborsOf(v, [&](node w, node u) {
				if (G.hasEdge(w, u)) {
					n++;
				}
			});

			mtx.lock();
			numerator += 2 * n;
			denominator += d * (d - 1);
			mtx.unlock();
		}
	});

	return (double) numerator / denominator;
}

} /* namespace NetworKit */