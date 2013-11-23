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


  std::pair<count,count> ecc(std::vector<node> array,int size) // returns the maximum entry of an unsorted array and its index
  {
    count j;
    count distance_max=0;
    for(int i=0; i<= size-1;i++)
    {
      if(array[i] > distance_max)
        {
          distance_max =array[i];
          j=i;
        }
    }
    return std::make_pair (distance_max,j);

  }

std::pair<count, count> GraphProperties::estimateDiameter_ck(const Graph& G) {
  count lowerBound = 0;
  count upperBound = std::numeric_limits<count>::max();

  int n = G.numberOfNodes();
  int maxi = minMaxDegree(G).second;
  int i = 1;
  int j = 0;
  std::vector<node> distances;
  std::vector<node> distance;
  
  std::vector<int> degree[maxi+1];
  count max;
  count v;
  
  G.forNodes([&](node u)
    {
      degree[G.degree(u)].push_back(u);           
    });
  int k = degree[maxi].back();
  std::cout<<k<<std::endl;

  while((upperBound- lowerBound) > 5)
  {
    while((degree[i].empty())&&(i <= (maxi/2)+1))               // improving lowerbound by computing ecc for a node with smallest degree       
    {
      i++;
    }
    node u = degree[i].back();
    distances = BFS().run(G,u);
    count ecc_result = ecc(distances,distances.size()).first; // yields ecc(u)
    if(ecc_result >= lowerBound)
      lowerBound = ecc_result;
    degree[i].pop_back();
    std::cout<< degree[i].size()<<std::endl;
    std::cout<<"lowerBound="<<lowerBound<<std::endl;

    while(degree[maxi-j].empty())          //improving upperbound by computing the diameter of spanning                           tree with root = node with highest degree.
    {
      j++;
    }
    distances = BFS().run(G,degree[maxi-j].back());
    max = ecc(distances,n).first;
    if(max >= lowerBound)
      lowerBound = max;
    std::cout<<"max="<< max <<std::endl;
    v = ecc(distances,n).second;
    std::cout<<"v="<< v <<std:: endl;
    distance = BFS().run(G,v);
    if((ecc(distance,n).first)<= upperBound)
      upperBound = ecc(distance,n).first;
    std::cout<<"upperBound="<< upperBound<<std::endl;
    degree[maxi-j].pop_back();
  }
  return std::make_pair(lowerBound, upperBound);
}

} /* namespace NetworKit */
