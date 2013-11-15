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
ClusteringCoefficient::approxAvgLocal(Graph& G) const
{
  // TODO implement approximative algorithm for average local cc
	return 0.0;
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
