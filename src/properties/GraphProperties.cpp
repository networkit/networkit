/*
 * GraphProperties.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "GraphProperties.h"
#include "../graph/BFS.h"
#include "../graph/BFSTree.h"
#include <random>
#include <algorithm>


namespace NetworKit {

GraphProperties::GraphProperties() {

}

GraphProperties::~GraphProperties() {

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

double GraphProperties::approximateGlobalClusteringCoefficient(
		const Graph& G, double approximationError, double probabilityError) {
	ApproximateClusteringCoefficient_Hoske acc;
	count numIters = acc.niters(approximationError, probabilityError);
	return acc.calculate(true, G, numIters);
}

std::pair<count, count> GraphProperties::estimatedDiameterRange_Feist(const Graph& G, count p) {
    count lowerBound = 0;
    count upperBound = std::numeric_limits<count>::max();

    node startNode;
    count inc = 0;
    count infDist = std::numeric_limits<count>::max();
    count n = G.numberOfNodes();       
    std::pair<std::vector<count>, node> resultOfBFS;

    std::default_random_engine rand;
    rand.seed(std::random_device()());
    std::uniform_int_distribution<node> range(0, n-1);

    std::vector<count> highToLow(n);
    std::iota(begin(highToLow), end(highToLow), 0);
    std::sort(begin(highToLow), end(highToLow), [&] (node v1, node v2) {
            return G.degree(v1) > G.degree(v2);
        });
  
    // check if G is a connected Graph
    {
        startNode = range(rand);

        BFS bfs;    
        resultOfBFS = bfs.run_Feist(G, startNode); 

        for(count& e : resultOfBFS.first){
            if (e == infDist) {
                lowerBound = infDist;
                upperBound = infDist;                    
                return std::make_pair(lowerBound, upperBound);
            }
        }        
    }

    while (upperBound - lowerBound >= p) {

        // 1. compute double sweep lower bound
        {
            startNode = range(rand);

            BFS bfs;    
            resultOfBFS = bfs.run_Feist(G, startNode); 
            node max_distance_node = resultOfBFS.second;             
                 
            BFS bfs2;
            resultOfBFS = bfs.run_Feist(G, max_distance_node);
                        
            if (resultOfBFS.first[resultOfBFS.second] > lowerBound)
                lowerBound = resultOfBFS.first[resultOfBFS.second];
    
        } // (1)

        // 2. compute tree upper bound 
        {        
            //startNode = d(e);
            startNode = highToLow[inc];
            std::vector<count> distances(n, infDist);
            std::queue<node> q;     
   
            distances[startNode] = 0;
            q.push(startNode);
        
            Graph spanningTree(n);

            while (! q.empty()) {
                node current = q.front();
                q.pop();
                        
                G.forNeighborsOf(current, [&](node neighbor) {
                        if (distances[neighbor] == infDist) {
                            q.push(neighbor);
                            distances[neighbor] = distances[current] + 1;
                            spanningTree.addEdge(current, neighbor);
                            startNode = neighbor;
                        }                               
                    });
            }                             
            
            BFS bfs;    
            resultOfBFS = bfs.run_Feist(spanningTree, startNode);    

            if (resultOfBFS.first[resultOfBFS.second] < upperBound) 
                upperBound = resultOfBFS.first[resultOfBFS.second];          
            
            inc++;
    
        } // (2)

    } //while
   
    return std::make_pair(lowerBound, upperBound);
}

count GraphProperties::DiameterRange_Feist(const Graph& G) {
    
    count diameter = 0;
    count current = 0;

    GraphDistance dist;

    G.forNodePairs([&](node u, node v){
            current = dist.unweightedDistance(G, u, v);
            if (current > diameter)
                diameter = current;
	});    


    return diameter;
}

std::pair<node, count> GraphProperties::eccentricity_Hoske(const Graph& G, node u) {
    static BFS bfs;
    auto dists = bfs.run(G, u);
    auto max_iter = max_element(begin(dists), end(dists));
    return {distance(begin(dists), max_iter), *max_iter};
}

std::pair<count, count> GraphProperties::estimatedDiameterRange_Hoske(const Graph& G, double error) {
	using namespace std;

	/* BFS that calls f with the visited edges and returns the node with largest distance from u. */
	/* Note: the function Graph::breadthFirstEdgesFrom that should
	   do the same has not been implemented! */
	auto bfs_edges = [&] (const Graph& G, node u, function<void(node, node)> f) -> node {
		queue<node> q;
		vector<bool> visited(G.numberOfNodes(), false);
		q.push(u);
		visited[u] = true;

		node x = u;
		while (!q.empty()) {
			x = q.front(); q.pop();
			G.forNeighborsOf(x, [&] (node y) {
				if (!visited[y]) {
					f(x, y);
					visited[y] = true;
					q.push(y);
				}
			});
		}
		return x;
	};

	/* Diameter estimate: lowerBounds <= diam(G) <= upperBound. */
	count lowerBound = 0;
	count upperBound = numeric_limits<count>::max();
	const count n = G.numberOfNodes();

	/* Nodes sorted decreasingly by degree. */
	vector<node> high_deg(n);
	iota(begin(high_deg), end(high_deg), 0);
	sort(begin(high_deg), end(high_deg), [&] (node u, node v) {
		return G.degree(u) > G.degree(v);
	});

	/* Random node. */
	static const default_random_engine random;
	auto random_node = bind(uniform_int_distribution<node>(0, n - 1), random);

	/* While not converged: update estimate. */
	count niter = 0;
	while ((upperBound - lowerBound) >= error*lowerBound && niter < n) {
		count ecc;

		/* ecc(u) <= diam(G) */
		node u = random_node();
		tie(ignore, ecc) = eccentricity_Hoske(G, u);
		lowerBound = max(lowerBound, ecc);

		/* diam(G) <= diam(BFS_Tree(v)) */
		node v = high_deg[niter];
		Graph bfs_tree(n);
		node w = bfs_edges(G, v, [&] (node a, node b) {
			bfs_tree.addEdge(a, b);
		});
		lowerBound = max(lowerBound, ecc);
		/* diam(T) = ecc_T(w) by problem 4. */
		tie(ignore, ecc) = eccentricity_Hoske(bfs_tree, w);
		upperBound = min(upperBound, ecc);

		niter++;
	}

	return {lowerBound, upperBound};
}


// returns the maximum entry of an unsorted array and its index
static std::pair<count, count> ecc(std::vector<node> distances) {
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

std::pair<count, count> GraphProperties::estimatedDiameterRange_Brueckner(
		const Graph& G) {

    using namespace std;

    count infDist = numeric_limits<count>::max();

    default_random_engine e;
    e.seed(random_device()());

    // Pick nodes randomly, weighted by their degree.
    vector<count> nodeWeights(G.numberOfNodes());
    G.forNodes([&](node v) {
            nodeWeights[v] = G.degree(v);
        });

	count lowerBound = 0;
	count upperBound = infDist;

    for (count i = 0;
         i < G.numberOfNodes() && upperBound - lowerBound > 5;
         i++) {

        discrete_distribution<node> nodeDistribution(nodeWeights.begin(), nodeWeights.end());
        node v = nodeDistribution(e);
        nodeWeights[v] = 0;
        BFSTree T(G, v);
        if (!T.spanning())
            return make_pair(infDist, infDist);
        else {
            lowerBound = max(lowerBound, T.depth());
            upperBound = min(upperBound, 2 * T.depth());
            BFS bfs;
            std::vector<count> dists = bfs.run(T, T.deepest());
            count maxDist = 0;
            for (auto dist : dists)
                maxDist = max(maxDist, dist);
            upperBound = min(upperBound, maxDist);
        }
    }

	return make_pair(lowerBound, upperBound);
}

count GraphProperties::exactDiameter_Brueckner(
        const Graph& G) {

    using namespace std;

    count diameter = 0;

	G.forNodesInRandomOrder([&](node v){
            BFS bfs;
            vector<count> distances = bfs.run(G, v);
            for (auto distance : distances) {
                if (diameter < distance)
                    diameter = distance;
            }
    });

    return diameter;
}

} /* namespace NetworKit */
