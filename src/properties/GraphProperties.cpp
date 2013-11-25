/*
 * GraphProperties.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "GraphProperties.h"


namespace NetworKit {

GraphProperties::GraphProperties() {
	// TODO Auto-generated constructor stub

}

GraphProperties::~GraphProperties() {
	// TODO Auto-generated destructor stub
}

std::vector<count> GraphProperties::degreeDistribution(const Graph& G) {
	count maxDegree = minMaxDegree(G).second;
	std::vector<count> distribution(maxDegree+1, 0);
	G.forNodes([&](node v){
		count i = G.degree(v);
		distribution[i]++;
	});
	return distribution;
}


std::vector<double> GraphProperties::localClusteringCoefficients(const Graph& G) {
	count n = G.numberOfNodes();
	std::vector<double> numerator(n); //
	std::vector<double> denominator(n); // $\deg(u) \cdot ( \deg(u) - 1 )$
	std::vector<double> coefficient(n); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

	G.parallelForNodes([&](node u){
		count edgeCount = 0;
		G.forEdgesOf(u, [&](node u, node v) {
			G.forEdgesOf(v, [&](node v, node w){
				if (G.hasEdge(u, w)) {
					edgeCount += 1;
				}
			});
		});

		numerator[u] = edgeCount; // factor 2 is omitted because each edge has been counted twice
	});

	G.parallelForNodes([&](node u){
		denominator[u] = G.degree(u) * (G.degree(u) - 1);
	});

	G.parallelForNodes([&](node u){
		if (denominator[u] == 0.0) {
			coefficient[u] = 0.0;
		} else {
			coefficient[u] = numerator[u] / denominator[u];
		}
	});

	return coefficient;
}

std::vector<double> GraphProperties::localClusteringCoefficientPerDegree(const Graph& G) {

	std::vector<count> degDist = degreeDistribution(G);
	std::vector<double> coefficient;
	std::vector<double> perDegree(degDist.size(), 0.0);

	if (G.numberOfNodes() > 1 ) {
		coefficient = localClusteringCoefficients(G);

			G.forNodes([&](node u){
				perDegree[G.degree(u)] += coefficient[u];
			});

			// get the average local clustering coefficient for nodes of each degreee
			for (index i = 2; i < degDist.size(); ++i) {
				if (degDist[i] == 0) {
					perDegree[i] = 0.0; // TODO: should this be -1
				} else {
					perDegree[i] = perDegree[i] / (double) degDist[i];
				}
			}
	}



	// allows to avoid the situation, when local clustering coefficient is calculated for 0-1 degree nodes.
	// These nodes are warranted not to be triangle centers, thus we avoid calculating the coefficients for the,
	degDist[0] = 0;
	if (G.numberOfNodes() > 0 ) degDist[1] = 0;

	return perDegree;
}

double GraphProperties::averageLocalClusteringCoefficient(const Graph& G) {
	std::vector<double> coefficients = GraphProperties::localClusteringCoefficients(G);
	double sum = 0.0;
	for (double c : coefficients) {
		sum += c;
	}
	double avg = sum / G.numberOfNodes();
	return avg;
}

std::pair<count, count> GraphProperties::minMaxDegree(const Graph& G) {

	count min = G.numberOfNodes();
	count max = 0;

	G.forNodes([&](node v){
		count d = G.degree(v);
		if (d < min) {
			min = d;
		}
		if (d > max) {
			max = d;
		}
	});

	return std::pair<count, count>(min, max);
}


double GraphProperties::averageDegree(const Graph& G) {
	count n = G.numberOfNodes();

	count degSum = G.parallelSumForNodes([&](node v){
		return G.degree(v);
	});

	double avgDeg = degSum / (double) n;
	return avgDeg;
}

// returns the maximum entry of an unsorted array and its index
std::pair<count, count> ecc(std::vector<node> distances) {
  count max_distance = std::numeric_limits<count>::min();
  node max_distance_node = 0;

  for(int i = 0; i < distances.size(); i++) {
    if(distances[i] > max_distance) {
      max_distance = distances[i];
      max_distance_node = i;
    }
  }
  return std::make_pair(max_distance, max_distance_node);
}

std::pair<count, count> GraphProperties::estimateDiameter_ck(const Graph& G) {
  count lowerBound = 0;
  count upperBound = std::numeric_limits<count>::max();

  int n = G.numberOfNodes();
  int maxDegree = minMaxDegree(G).second;

  std::vector<node> nodesWithDegree[maxDegree + 1];
  
  G.forNodes([&](node u) {
    nodesWithDegree[G.degree(u)].push_back(u);           
  });

  int i = 1, j = 0;
  while(upperBound - lowerBound > 5) {
    // improving lower bound by computing ecc for a node with smallest degree
    while(nodesWithDegree[i].empty() && i <= maxDegree/2 + 1) {
      i++;
    }
    node u = nodesWithDegree[i].back();
    std::vector<node> distances = BFS().run(G, u);
    count ecc_result = ecc(distances).first; // yields ecc(u)
    if(ecc_result > lowerBound) {
      lowerBound = ecc_result;
    }
    nodesWithDegree[i].pop_back();

    // improving upper bound by computing the diameter of spanning
    // tree with root = node with highest degree.
    while(nodesWithDegree[maxDegree - j].empty()) {
      j++;
    }
    u = nodesWithDegree[maxDegree - j].back();
    distances = BFS().run(G, u);
    auto ecc_pair = ecc(distances);
    ecc_result = ecc_pair.first;
    if(ecc_result > lowerBound) {
      lowerBound = ecc_result;
    }

    Graph T(n);
    G.forNodes([&](node w) {
      G.forEdgesOf(w, [&](node u, node v) {
        if(distances[u] == distances[v] + 1) {
          T.addEdge(u,v);
        }
      });
    });
    u = ecc_pair.second;
    distances = BFS().run(T, u);
    ecc_result = ecc(distances).first;
    if(ecc_result < upperBound) {
      upperBound = ecc_result;
    }
    nodesWithDegree[maxDegree - j].pop_back();
  }
  return std::make_pair(lowerBound, upperBound);
}

} /* namespace NetworKit */
