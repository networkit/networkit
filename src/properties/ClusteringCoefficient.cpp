/*
 * ClusteringCoefficient.cpp
 *
 *  Created on: 08.04.2013
 *      Author: cls
 */

#include "ClusteringCoefficient.h"

namespace NetworKit {

ClusteringCoefficient::ClusteringCoefficient()
{
	// TODO Auto-generated constructor stub

}

ClusteringCoefficient::~ClusteringCoefficient()
{
	// TODO Auto-generated destructor stub
}

double
ClusteringCoefficient::avgLocal(Graph& G) const
{

	count n = G.numberOfNodes();
	std::vector<double> coefficient(n); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$


	G.parallelForNodes([&](node u){
		count triangles = 0;
		G.forEdgesOf(u, [&](node u, node v) {
			G.forEdgesOf(v, [&](node v, node w){
				if (G.hasEdge(u, w)) {
					triangles += 1;
				}
			});
		});

    count d = G.degree(u);
    coefficient[u] = (double)triangles / (double)(d * (d - 1));
	});

	double cc = G.parallelSumForNodes([&](node u){
		return coefficient[u];
	});

	cc /= (double)n;

	return cc;
}

double
ClusteringCoefficient::approxAvgLocal(Graph& G, count tries) const
{
	count n = G.numberOfNodes();

	// WARNING: I assume RAND_MAX to be larger than n. If this should not hold for an application
	// or implementation of the standard library, a more sophisticated version of determining a 
	// vertex uniformly at random must be used.

	count triangles = 0;
	for (count k = 0; k < tries; ++k) {
		node v = rand() % n;

		if (G.degree(v) < 2) {
			// this iteration is ignored, since this vertex can never be part of a triangle,
			// nor middle point of a path of length 3

			--k;
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

	return triangles / tries;
}


double
ClusteringCoefficient::global(Graph& G) const
{

	count n = G.numberOfNodes();
	std::vector<count> triangles(n); // triangles including node u (every triangle is counted six times)
	std::vector<count> triples(n); // triples around node u


	G.parallelForNodes([&](node u){
		count trias = 0;
		G.forEdgesOf(u, [&](node u, node v) {
			G.forEdgesOf(v, [&](node v, node w){
				if (G.hasEdge(u, w)) {
					trias += 1;
				}
			});
		});
    
    triangles[u] = tris;
	});
  
  double denominator = G.parallelSumForNodes([&](node u){
		return G.degree(u) * (G.degree(u) - 1);
	});

	double cc = G.parallelSumForNodes([&](node u){
		return triangles[u];
	});

	cc /= (2 * denominator); // factor 2 because triangles counted six times, but need them three times

	return cc;
}


double
ClusteringCoefficient::approxGlobal(Graph& G) const
{
  // TODO implement approximative algorithm for global cc
	return 0.0;
}

} /* namespace NetworKit */
