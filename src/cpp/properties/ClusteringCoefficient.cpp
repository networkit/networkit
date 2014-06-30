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

namespace NetworKit {

std::vector<double> ClusteringCoefficient::exactLocal(Graph &G) {
	count z = G.upperNodeIdBound();
	std::vector<double> coefficient(z); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

	G.balancedParallelForNodes([&](node u) {

#if 0
		// TODO:
		// for each vertex u
		// retrieve neighborhood of u
		// examine each pair (v, w) of neighbors of u
		// if (v, w) \in E then increase triangles

		count d = G.degree(u);
	    if (d < 2) {
	      coefficient[u] = 0.0;
	    } else {
		      count triangles = 0;
		      std::vector<node> neigh = G.neighbors(u);
		      for (index i = 0; i < d; ++i) {
		    	  for (index j = i+1; j < d; ++j) {
		    		  if (G.hasEdge(neigh[i], neigh[j])) {
		    			  triangles++;
		    		  }
		    	  }
		      }
		      coefficient[u] = (double) triangles / (double)(d * (d - 1)); // No division by 2 since triangles are counted twice as well!
	    }
#endif


		count d = G.degree(u);
		std::unordered_set<node> uNeighbors; // set for O(1) time access to u's neighbors
		G.forNeighborsOf(u, [&](node v){
			if (v != u) {
				uNeighbors.insert(v);
			}
		});
	    
	    if (d < 2) {
	      coefficient[u] = 0.0;
	    } else {
	      count triangles = 0;
	      G.forEdgesOf(u, [&](node u, node v) {
	        G.forEdgesOf(v, [&](node v, node w){
	          if (uNeighbors.find(w) != uNeighbors.end()) { // w is also a neighbor of u
	            triangles += 1;
	          }
	        });
	      });
	      coefficient[u] = (double) triangles / (double)(d * (d - 1)); // No division by 2 since triangles are counted twice as well!
	    }
	});

	return coefficient;
}

double ClusteringCoefficient::avgLocal(Graph& G) {

	auto coefficients = exactLocal(G); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

	double sum = 0.0;
	count size = 0;

	G.forNodes([&](node u) {
		if (G.degree(u) >= 2) {
			sum += coefficients[u];
			size++;
		}
	});

	return sum / (double) size;
}

double ClusteringCoefficient::approxAvgLocal(Graph& G, const count trials) {

	double triangles = 0;
	for (count k = 0; k < trials; ++k) {
		node v = G.randomNode();
		TRACE("trial ", k, " sampled node ", v);

		if (G.degree(v) < 2) {
			// this vertex can never be part of a triangle,
			// nor middle point of a path of length 3
			--k;  // do not count trial
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
	std::vector<count> triples(z); // triples around node u


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
