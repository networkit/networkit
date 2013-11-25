/*
 * GraphProperties.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "GraphProperties.h"
#include "random"


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


} /* namespace NetworKit */

