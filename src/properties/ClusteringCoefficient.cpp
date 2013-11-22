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

std::vector<double>
ClusteringCoefficient::local(Graph &G) const
{
	count n = G.numberOfNodes();
	std::vector<double> coefficient(n); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

	G.parallelForNodes([&](node u){
	    count d = G.degree(u);
	    
	    if (d < 2) {
	      coefficient[u] = 0.0;
	    } else {
	      count triangles = 0;
	      G.forEdgesOf(u, [&](node u, node v) {
	        G.forEdgesOf(v, [&](node v, node w){
	          if (G.hasEdge(u, w)) {
	            triangles += 1;
	          }
	        });
	      });
	      coefficient[u] = (double)triangles / (double)(d * (d - 1));
	    }
	});

	return coefficient;
}

double
ClusteringCoefficient::avgLocal(Graph& G) const
{
	count n = G.numberOfNodes();
	std::vector<double> coefficient(n); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

	coefficient = this->local(G);

	double cc = G.parallelSumForNodes([&](node u){
		return coefficient[u];
	});

	cc /= (double)n;

	return cc;
}

double
ClusteringCoefficient::approxAvgLocal(Graph& G, const count tries) const
{
	count n = G.numberOfNodes();

	// WARNING: I assume RAND_MAX to be larger than n. If this should not hold for an application
	// or implementation of the standard library, a more sophisticated version of determining a 
	// vertex uniformly at random must be used.

	double triangles = 0;
	for (count k = 0; k < tries; ++k) {
		node v = rand() % n;

		if (G.degree(v) < 2) {
			// this vertex can never be part of a triangle,
			// nor middle point of a path of length 3
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

	return triangles / (double)tries;
}


double
ClusteringCoefficient::global(Graph& G) const
{
	count n = G.numberOfNodes();
	std::vector<count> triangles(n); // triangles including node u (every triangle is counted six times)
	std::vector<count> triples(n); // triples around node u


	G.parallelForNodes([&](node u){
		count tr = 0;
    	if (G.degree(u) > 1) {
     	 	G.forEdgesOf(u, [&](node u, node v) {
       			G.forEdgesOf(v, [&](node v, node w){
         			if (G.hasEdge(u, w)) {
            			tr += 1;
          			}
        		});
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


double
ClusteringCoefficient::approxGlobal(Graph& G, const count tries) const
{
	count n = G.numberOfNodes();

  // Calculate prefix sum over the nodes where each node v counts deg(v)*(deg(v)-1) times
	std::vector<count> weight(n);
	count psum = 0;
	for (index i = 0; i < n; i++) {
		psum += G.degree(i) * (G.degree(i) - 1);
		weight[i] = psum;
	}

	// WARNING: I assume RAND_MAX to be larger than PSUM. If this should not hold for an application
	// or implementation of the standard library, a more sophisticated version of determining a 
	// vertex uniformly at random must be used.

	double triangles = 0;
	for (count k = 0; k < tries; ++k) {
		count r = rand() % (psum+1);

		// plain old binary search:
		index low = 0; 
		index high = G.numberOfNodes();
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

		if (G.degree(v) < 2) {
			// this vertex can never be part of a triangle,
			// nor middle point of a path of length 3
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

} /* namespace NetworKit */
