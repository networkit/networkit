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

	count numerator = 0;
	count denominator = 0;

	G.parallelForNodes([&](node v) {
		const count d = G.degree(v);
		if (d >= 2) {
			count n = 0;

			// count triangles of v
			G.forNeighborsOf(v, [&](node w) {
				G.forNeighborsOf(v, [&](node u) {
					if (w != u && G.hasEdge(w, u)) {
						n++; // this will count every triangle twice because we have the pair (w,u) and (u,w)
					}
				});
			});

			mtx.lock();
			numerator += n; // no times 2 because we counted every triangle twice
			denominator += d * (d - 1);
			mtx.unlock();
		}
	});

	return (double) numerator / denominator;
}

} /* namespace NetworKit */