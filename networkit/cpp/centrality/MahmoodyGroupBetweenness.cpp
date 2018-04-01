/*
 * MahmoodyGroupBetweenness.h
 *
 *  Created on: 13.03.2018
 *      Author: Marvin Pogoda
 */

#include "MahmoodyGroupBetweenness.h"
#include <omp.h>
#include "../distance/BFS.h"
#include "../distance/SSSP.h"
#include "../auxiliary/BucketPQ.h"
#include <math.h>
#include "../auxiliary/Random.h"  


namespace NetworKit{

MahmoodyGroupBetweenness::MahmoodyGroupBetweenness(const Graph& G,count groupSize,double epsilon): G(G), groupSize(groupSize),epsilon(epsilon)
{
}


void MahmoodyGroupBetweenness::run()
{	
	
	
	//Create data structures for the hypergraph.
	std::vector<int64_t> v(G.numberOfNodes());
	std::vector<std::vector<node>> adjacencyList(G.numberOfNodes());
	std::vector<int> hyperEdges;	

	int samples = 2 * log(pow(G.numberOfNodes(),3)) / epsilon; 
	Aux::BucketPQ nodeDegrees(v,-samples,1);
	omp_lock_t lock;
	omp_init_lock(&lock);
	#pragma omp parallel for
	for(int l = 0;l < samples;l++)   
	{
		node s= G.randomNode();
		node t;
		do{
			t= G.randomNode();
			
		}while(s == t);

		BFS bfs(G,s,true,true,t);
		bfs.run();
		std::set<std::vector<node>> shortestPaths = bfs.getPaths(t);

		if(shortestPaths.size() == 0){		//If the selected nodes are in different connected components, the hyperedge is an empty set.
		continue;				//Chooseing nodes in different connected components wont affect the algorithm.
		}					//(See Mahmoody "Scalable Betweenness Centrality Maximization via Sampling",page 4,Lemma 3,2016)
			
		//Uniformly select a shortest path
		std::set<std::vector<node>>::const_iterator iterator(shortestPaths.begin());
		advance(iterator,Aux::Random::integer(shortestPaths.size() - 1));
		std::vector<node> newHyperEdge = *iterator;
		omp_set_lock(&lock);
		
		//Insert length of hyperedge and save the position for the ajdacency list.
		hyperEdges.push_back(newHyperEdge.size());
		node hyperEdgeStart= hyperEdges.size() - 1;
		
		for(node const& n:newHyperEdge){
			//Updates the degree of nodes in the hypergraph
			nodeDegrees.changeKey(nodeDegrees.getKey(n) - 1,n);
			
			//Insert new Hyperedge into the hyperedge-list
			hyperEdges.push_back(n);

			//Update AdjancecyList
			adjacencyList[n].push_back(hyperEdgeStart);
		}
		omp_unset_lock(&lock);
	}
	
    	//Extract nodes with highest degress.
	for(count j = 0;j < groupSize;j++){
	std::pair<int,node> elem = nodeDegrees.extractMin();
	node n = elem.second;
	
	//Update degrees
	for(auto hyperEdge:adjacencyList[n]){
		count start= hyperEdge+1;
		count end = start+hyperEdges[hyperEdge];
		for(count i=start;i < end;i++){
			if(hyperEdges[i] != n){
				node decreasedNode = hyperEdges[i];
				nodeDegrees.changeKey(nodeDegrees.getKey(decreasedNode)+1,decreasedNode);
			}
		}
		hyperEdges[start - 1] = 0;
	}
	maxGroup.push_back(n);
	}
	
	hasRun=true;
}
	
}/* namespace NetworKit */

