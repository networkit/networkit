/*
 * ClusteringCoefficient.cpp
 *
 *  Created on: 08.04.2013
 *      Author: Lukas Barth, David Weiss
 */

#include <unordered_set>
 
#include "ClusteringCoefficient.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/Log.h"
#include <omp.h>

namespace NetworKit {

std::vector<double> ClusteringCoefficient::exactLocal(Graph &G) {
	count z = G.upperNodeIdBound();
	std::vector<double> coefficient(z); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

	std::vector<std::vector<bool> > nodeMarker(omp_get_max_threads());

	for (auto & nm : nodeMarker) {
		nm.resize(z, false);
	}

	G.balancedParallelForNodes([&](node u) {
		count d = G.degree(u);

		if (d < 2) {
			coefficient[u] = 0.0;
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

			coefficient[u] = (double) triangles / (double)(d * (d - 1)); // No division by 2 since triangles are counted twice as well!
		}
	});

	return coefficient;
}

double ClusteringCoefficient::avgLocal(Graph& G) {

	auto coefficients = exactLocal(G); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

	double sum = 0.0;

	G.forNodes([&](node u) {
		sum += coefficients[u];
	});

	return sum / (double) G.numberOfNodes();
}

double ClusteringCoefficient::approxAvgLocal(Graph& G, const count trials) {

	double triangles = 0;
	for (count k = 0; k < trials; ++k) {
		node v = G.randomNode();
		TRACE("trial ", k, " sampled node ", v);

		if (G.degree(v) < 2) {
			// this vertex can never be part of a triangle,
			// (nor middle point of a path of length 3)
			continue;
		}

		TRACE("deg(v) = ", G.degree(v));
		node u = G.randomNeighbor(v);
		node w = G.randomNeighbor(v);
		TRACE("u=", u);
		TRACE("w=", w);

		// TODO This could be sped up for degree(v) == 2...
		while (u == w) {
			w = G.randomNeighbor(v);
			TRACE("w=", w);
		}

		if (G.hasEdge(u,w)) {
			triangles++;
		}
	}

	return triangles / (double) trials;
}


double ClusteringCoefficient::exactGlobal(Graph& G) {
	count z = G.upperNodeIdBound();
	std::vector<count> triangles(z); // triangles including node u (every triangle is counted six times)

	std::vector<std::vector<bool> > nodeMarker(omp_get_max_threads());
	for (auto &nm : nodeMarker) {
		nm.resize(z, false);
	}

	G.balancedParallelForNodes([&](node u){

		size_t tid = omp_get_thread_num();
		count tr = 0;

		if (G.degree(u) > 1) {
			G.forEdgesOf(u, [&](node u, node v) {
				nodeMarker[tid][v] = true;
			});

			G.forEdgesOf(u, [&](node u, node v) {
				G.forEdgesOf(v, [&](node v, node w) {
					if (nodeMarker[tid][w]) {
						tr += 1;
					}
				});
			});

			G.forEdgesOf(u, [&](node u, node v) {
				nodeMarker[tid][v] = false;
			});
		}

		triangles[u] = tr;
	});
  
  double denominator = G.parallelSumForNodes([&](node u){
		return G.degree(u) * (G.degree(u) - 1);
	});

	double cc = G.parallelSumForNodes([&](node u){
		return triangles[u];
	});

	cc /= denominator; 

	return cc;
}


double ClusteringCoefficient::approxGlobal(Graph& G, const count trials) {
	count z = G.upperNodeIdBound();

  // Calculate prefix sum over the nodes where each node v counts deg(v)*(deg(v)-1) times
	std::vector<count> weight(z);
	count psum = 0;
	G.forNodes([&](node v) {
		psum += G.degree(v) * (G.degree(v) - 1);
		weight[v] = psum;
	});

	// WARNING: I assume RAND_MAX to be larger than PSUM. If this should not hold for an application
	// or implementation of the standard library, a more sophisticated version of determining a 
	// vertex uniformly at random must be used.

	double triangles = 0;
	for (count k = 0; k < trials; ++k) {
		count r = Aux::Random::integer(psum - 1);

		// plain old binary search:
		index low = 0; 
		index high = G.upperNodeIdBound();
		while (low < high) {
			index middle = (low + high) / 2;

			if (weight[middle] < r) {
				low = middle + 1;
			} else if (weight[middle] > r) {
				high = middle;
			} else {
				low = high = middle;
			}
		}

		node v = low; // ewww.. setting a vertex to an index.. but works.

		count vDeg = G.degree(v);
		if (vDeg < 2) {
			// this vertex can never be part of a triangle,
			// nor middle point of a path of length 3
			--k; // do not count trial
			continue;
		}

		node u = G.randomNeighbor(v);
		node w = G.randomNeighbor(v);

		// TODO This could be sped up for degree(v) == 2...
		while (u == w) {
			w = G.randomNeighbor(v);
		}

		if (G.hasEdge(u,w)) {
			triangles++;
		}
	}

	return triangles / trials;
}

} /* namespace NetworKit */
