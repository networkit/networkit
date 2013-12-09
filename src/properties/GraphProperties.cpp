/*
 * GraphProperties.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include <list>

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

	std::pair<count,count> ecc(std::vector<node> array,int size) 
	{
		count infDist = std::numeric_limits<count>::max();
		count distance_max=0;
		count j=0;
		for(int i=0; i<= size-1;i++)
		{
			if((array[i] > distance_max)&&(array[i]<infDist))
				{
				  distance_max = array[i];
				  j=i;
				}
		}
		return std::make_pair(distance_max,j);

	}

std::pair<std::vector<double>,double> normed_OckerReichard(std::vector<double> nextVector, int n){
double norm=0;
for(int i=0; i<= n-1; i++)
	{
		norm = norm + (nextVector[i]*nextVector[i]);
	}
	norm = (double) sqrt(norm);
	//std::cout<<"norm="<<norm<<std::endl;
	for(int i=0; i<= n-1; i++)
	{
		nextVector[i]=nextVector[i]/norm;
		
	}
	return make_pair(nextVector,norm);
}

std::vector<double> Matrixmultiplication_OckerReichard(std::vector<std::vector<double>> adjacencymatrix,int n, std::vector<double> vec){

double sum=0;
std::vector<double> nextVector(n);
for(int i=0; i <= n-1; i++)
	{
		for(int j=0; j<= n-1; j++)
		{
			sum=sum+adjacencymatrix[i][j]*vec[j];
		}
	nextVector[i]=sum;
	sum=0;
	}
	return nextVector;

}

std::vector<double> GraphProperties::EVZ_OckerReichard(const Graph& G) {
int n= G.numberOfNodes();
std::vector<double> Vector_current(n);
std::vector<double> Vector_next(n);

std::vector<double> distance(n);
double distanc;
double epsilon= 0.001;
double sum=0;

for(int i=0; i<= n-1; i++)
{
	Vector_current[i]=1;
	Vector_next[i]=1;
}
std::cout<<"Push_backed"<<std::endl;
std::vector< std::vector<double> > adjacencymatrix(n, std::vector<double>(n));
for(node i=0; i<= n-1; i++)
{
	for(node j=0; j<= n-1; j++)
	{
		if(G.hasEdge(i,j))
			adjacencymatrix[i][j]=1;
		else{
			adjacencymatrix[i][j]=0;
		}
		
	}
}

do{
	Vector_current = Vector_next;
	Vector_next = Matrixmultiplication_OckerReichard(adjacencymatrix, n, Vector_current);
	Vector_next =normed_OckerReichard(Vector_next,n).first;
	for(int i=0; i<= n-1; i++)
	{
		//std::cout<<"Vector["<<i<<"]="<<Vector_next[i]<<std::endl;
	}
	for(int i=0; i<= n-1; i++)
	{
		distance[i]=Vector_next[i]-Vector_current[i];
		//std::cout<<"distance["<<i<<"]="<<distance[i]<<std::endl;
	}

	distanc=normed_OckerReichard(distance, n).second;
	//std::cout<<distanc<<std::endl;
  }while(distanc>epsilon);
return Vector_next;
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

std::vector<double> GraphProperties::betweennessCentrality_OckerReichard(const Graph& g) {
  const int n = g.numberOfNodes();
  std::vector<double> centrality(n, 0.0);

  g.forNodes([&](node s) {
    std::vector<std::vector<node>> predecessors(n, std::vector<node>());
    std::vector<int> distance(n, -1);
    std::vector<double> sigma(n, 0.0);
    std::vector<node> visitedNodes;
    std::list<node> nodesToVisit;

    distance[s] = 0;
    sigma[s] = 1.0;
    nodesToVisit.push_back(s);

    // BFS
    while(!nodesToVisit.empty()) {
      node v = nodesToVisit.front();
      nodesToVisit.pop_front();
      visitedNodes.push_back(v);

      g.forNeighborsOf(v, [&](node w) {
        if(distance[w] < 0) {
          // We are visiting w the first time
          nodesToVisit.push_back(w);
          distance[w] = distance[v] + 1;
        }

        if(distance[w] == distance[v] + 1) {
          sigma[w] += sigma[v];
          predecessors[w].push_back(v);
        }
      });
    }

    std::vector<double> delta(n, 0.0);

    // Traverse nodes in order of non-increasing distance
    while(!visitedNodes.empty()) {
      node w = visitedNodes.back();
      visitedNodes.pop_back();

      for(node v: predecessors[w]) {
        delta[v] += sigma[v] / sigma[w] * (1.0 + delta[w]);
      }

      if(w != s) {
        centrality[w] += delta[w];
      }
    }
  });

  return centrality;
}

} /* namespace NetworKit */
