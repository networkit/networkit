/*
 * SelectiveSCAN.cpp
 *
 *  Created on: 14.06.2013
 *      Author: cls
 */

#include "SelectiveSCAN.h"

namespace NetworKit {

SelectiveSCAN::SelectiveSCAN(): epsilon(0.7), mu(3) {

}

SelectiveSCAN::~SelectiveSCAN() {
	// TODO Auto-generated destructor stub
}

std::unordered_map<node, std::unordered_set<node>> SelectiveSCAN::run(Graph& G, std::unordered_set<node> set){

	std::unordered_map<node, node> nodesState;
	std::unordered_map<node, std::unordered_set<node>> communities;
	std::unordered_set<node> community;

	G.forNodes ([&](node u){
		nodesState.insert(std::pair<node,int>(u, -1));
	});
	for (node u : set) {
		std::pair<bool,std::vector<node>> isCore = this->isCore(u, G);
		std::vector<node> candidates;
		for(node v : isCore.second){
			candidates.push_back(v);
		}
		if ((nodesState.find(u))->second == -1  && isCore.first){
			expandCore(u, u, &community, &nodesState, &candidates, G);
		} else if ((nodesState.find(u))->second >= 0) {
			community = communities.find((nodesState.find(u))->second)->second;
		} else {
			bool clustered = false;
			while (!candidates.empty() && !clustered) {
				if(this->isCore(*candidates.begin(), G).first) {
					std::pair<bool,std::vector<node>> isCore = this->isCore(u, G);
					expandCore(*candidates.begin(), u, &community, &nodesState, &candidates, G);
					clustered = true;
				} else {
					candidates.erase(candidates.begin());
				}
			}
		}
		communities.insert({u, community});
	}
	return communities;
}

double SelectiveSCAN::nodeDistance(node u, node v, Graph& G) {

	int inter = 0;
	int uni = 0;
	G.forNeighborsOf(u, [&](node x){
		if (x != v ) {
			if (G.hasEdge(x, v)) {
				inter++;
				uni++;
			} else {
				uni++;
			}
		}
	});
	G.forNeighborsOf(v, [&](node x){
		if (x != u ) {
			if (!G.hasEdge(x, u)) {
				uni++;
			}
		}
	});
	return 1- ((double) (inter + 2)/ (double) (uni +2));
}

std::pair<bool,std::vector<node>> SelectiveSCAN::isCore(node u, Graph& G) {

	bool core = false;
	std::vector<node> similarNeighbors;
	int count = 0;
	G.forNeighborsOf(u, [&](node v){
		if (this->nodeDistance(u, v, G) <= this->epsilon) {
			count++;
			similarNeighbors.push_back(v);
		}
	});
	if (count >= this->mu) {
		core = true;
	}
	return std::pair<bool,std::vector<node>>(core, similarNeighbors);
}

void SelectiveSCAN::expandCore(node core, node label, std::unordered_set<node>* community,
		std::unordered_map<node, node>* nodesState, std::vector<node>* candidates, Graph& G) {

	std::pair<bool,std::vector<node>> isCore;
	community->insert(core);
	nodesState->find(core)->second = label;
	for (node v : *candidates) {
		nodesState->find(v)->second = label;
		community->insert(v);
		isCore = this->isCore(v, G);
		if (isCore.first) {
			for (node x : isCore.second) {
				int tmp = nodesState->find(x)->second;
				if (tmp < 0) {
					nodesState->find(x)->second = label;
					community->insert(x);
				}
				if (tmp == -1) {
					candidates->push_back(x);
				}
			}
		}
		candidates->erase(candidates->begin());
	}
}

} /* namespace NetworKit */
