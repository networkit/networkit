/*
 * BetweennessCentrality.cpp
 *
 *  Created on: 09.12.2013
 *      Author: Lukas Barth, David Weiß
 */

#include "BetweennessCentrality.h"

#include <iostream>

#include "../graph/BFS.h"

namespace GrauBart {
using namespace NetworKit;

BetweennessCentrality::BetweennessCentrality(Graph &G)
{
	this->g = &G;
}

BetweennessCentrality::~BetweennessCentrality()
{
	// TODO Auto-generated destructor stub
}

void 
BetweennessCentrality::computeTrees() 
{
	this->g->forNodes([&](node s){
		this->dists[s] = BFS().run(*(this->g), s);
	});
}

void
BetweennessCentrality::computeSuccessors()
{
	this->g->forNodes([&](node s){
		this->g->forNodes([&](node v){
			// Initialize data
			if (s != v) {
				this->numShortestPaths[s][v] = 0;
			} else {
				// Diagonal has to be set to 1, as there is 1 path v->v by default.
				this->numShortestPaths[s][v] = 1;
			}
		});
	});

	this->g->forNodes([&](node s){

		this->g->forEdges([&](node v, node w) {
			if (this->dists[s][v] == this->dists[s][w] + 1) {
				this->successors[s][w].push_front(v);
				this->numShortestPaths[s][v] += 1;
			} else if (this->dists[s][v] == this->dists[s][w] - 1) {
				this->successors[s][v].push_front(w);
				this->numShortestPaths[s][w] += 1;
			}
		});
	});
}

void
BetweennessCentrality::computeDependencies()
{
	this->g->forNodes([&](node s){
		// To start at the leaves, we'll have to sort ascending w.r.t. the BFS dist
		std::multimap<count, node> reverseDist;

		node v = 0; // this is sooooo dirty. I whished NetworKit would use maps instead of vectors :(
		for (auto it = this->dists[s].begin(); it != this->dists[s].end(); ++it) {
			reverseDist.insert(std::pair<count, node>(*it, v));
			v++; // ewwwwwww
		}

		// Now, iterate in reverse order to get the biggest distances first
		for (auto it = reverseDist.rbegin(); it != reverseDist.rend(); ++it) {
			node v = it->second;
			double dependency = 0.0;

			// For every vertex, iterate its successors
			for (auto succIt = this->successors[s][v].begin(); 
					succIt != this->successors[s][v].end(); ++succIt) {
				node w = *succIt;
				// we have to have settled the successor
				assert(this->dependencies[s].find(w) != this->dependencies[s].end());

				dependency += ((double)this->numShortestPaths[s][v] / (double)this->numShortestPaths[s][w]) 
							* (1 + this->dependencies[s][w]);
			}

			this->dependencies[s][v] = dependency;
		}
	});	
}

void
BetweennessCentrality::accumulateDependencies()
{	
	this->g->forNodes([&](node v){
		double BC = 0.0;
		this->g->forNodes([&](node s){
			BC += this->dependencies[s][v];
		});

		this->BCs[v] = BC;
	});
}

void
BetweennessCentrality::run() 
{
	this->computeTrees();
	this->computeSuccessors();
	this->computeDependencies();
	this->accumulateDependencies();
}

std::hash_map<node, double>
BetweennessCentrality::getCentrality() 
{
	return this->BCs;
}

} /* namespace GrauBart */
