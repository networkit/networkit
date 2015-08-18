#include "LocalClusteringCoefficient.h"
#include <omp.h>

namespace NetworKit {

LocalClusteringCoefficient::LocalClusteringCoefficient(const Graph& G) : Centrality(G, false, false) {
	if (G.isDirected()) throw std::runtime_error("Not implemented: Local clustering coefficient is currently not implemted for directed graphs");
	if (G.numberOfSelfLoops()) throw std::runtime_error("The Local Clustering Coeffiecient implementation cannot yet deal with self-loops");
}

void LocalClusteringCoefficient::run() {
	count z = G.upperNodeIdBound();
	scoreData.clear();
	scoreData.resize(z); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

	std::vector<std::vector<bool> > nodeMarker(omp_get_max_threads());

	for (auto & nm : nodeMarker) {
		nm.resize(z, false);
	}

	G.balancedParallelForNodes([&](node u) {
		count d = G.degree(u);

		if (d < 2) {
			scoreData[u] = 0.0;
		} else {
			size_t tid = omp_get_thread_num();
			count triangles = 0;

			G.forEdgesOf(u, [&](node u, node v) {
				nodeMarker[tid][v] = true;
			});

			G.forEdgesOf(u, [&](node u, node v) {
				G.forEdgesOf(v, [&](node v, node w) {
					if (nodeMarker[tid][w]) {
						triangles += 1;
					}
				});
			});

			G.forEdgesOf(u, [&](node u, node v) {
				nodeMarker[tid][v] = false;
			});

			scoreData[u] = (double) triangles / (double)(d * (d - 1)); // No division by 2 since triangles are counted twice as well!
		}
	});
	hasRun = true;
}


double LocalClusteringCoefficient::maximum() {
	return 1.0;
}

}
