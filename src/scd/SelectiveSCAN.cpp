/*
 * SelectiveSCAN.cpp
 *
 *  Created on: 14.06.2013
 *      Author: cls
 */

#include "SelectiveSCAN.h"

namespace NetworKit {

SelectiveSCAN::SelectiveSCAN(const Graph& G, NodeDistance& distMeasure, double epsilon, double mu): SelectiveCommunityDetector(G), distMeasure(&distMeasure), epsilon(epsilon), mu(mu) {
	DEBUG("preprocessing node distances");
	this->distMeasure->preprocess(); // distances depend on the graph, which is constant, and should be calculated only once

}

SelectiveSCAN::~SelectiveSCAN() {
	// TODO Auto-generated destructor stub
}

std::unordered_map<node, std::unordered_set<node>> SelectiveSCAN::run(std::unordered_set<node> set){

	std::unordered_map<node, node> nodesState;
	std::unordered_map<node, std::unordered_set<node>> communities;
	std::unordered_set<node> community;

	G.forNodes ([&](node u){
		nodesState.insert(std::pair<node,int>(u, -1));
	});
	for (node u : set) {
		std::pair<bool,std::unordered_set<node>> isCore = this->isCore(u);
		if ((nodesState.find(u))->second == -1  && isCore.first){
			expandCore(u, u, &community, &nodesState, &isCore.second);
		} else if (((nodesState.find(u))->second != -1) && ((nodesState.find(u))->second != -2)) {
			community = communities.find((nodesState.find(u))->second)->second;
		} else {
			bool clustered = false;
			nodesState.find(u)->second = -2;
			std::pair<bool,std::unordered_set<node>> isCore = this->isCore(u);
			std::pair<bool,std::unordered_set<node>> candidates;

			while (!isCore.second.empty() && !clustered) {
				node core = *isCore.second.begin();
				candidates = this->isCore(core);
				if(candidates.first) {
					expandCore(core, u, &community, &nodesState, &candidates.second);
					clustered = true;
				} else {
					isCore.second.erase(isCore.second.begin());
				}
			}
		}
		communities.insert({u, community});
		community.clear();
	}
	return communities;
}

//double SelectiveSCAN::nodeDistance(node u, node v) {
//
//	int inter = 0;
//	int uni = 0;
//	G.forNeighborsOf(u, [&](node x){
//		if (x != v && x!= u ) {
//			if (G.hasEdge(x, v)) {
//				inter++;
//				uni++;
//			} else {
//				uni++;
//			}
//		}
//	});
//	G.forNeighborsOf(v, [&](node x){
//		if (x != u && x != v) {
//			if (!G.hasEdge(x, u)) {
//				uni++;
//			}
//		}
//	});
//	return 1- ((double) (inter + 2)/ (double) (uni +2));
//}

std::pair<bool,std::unordered_set<node>> SelectiveSCAN::isCore(node u) {
	bool core = false;
	std::unordered_set<node> similarNeighbors;
	int count = 0;
	G.forNeighborsOf(u, [&](node v){
		if (this->distMeasure->distance(u, v) <= this->epsilon) {
			count++;
			similarNeighbors.insert(v);
		}
	});
	if (count >= this->mu) {
		core = true;
	}

	return std::pair<bool,std::unordered_set<node>>(core, similarNeighbors);
}

void SelectiveSCAN::expandCore(node core, node label, std::unordered_set<node>* community,
		std::unordered_map<node, node>* nodesState, std::unordered_set<node>* candidates) {
	std::pair<bool,std::unordered_set<node>> isCore;
	community->insert(core);
	nodesState->find(core)->second = label;

	while (!candidates->empty()) {
		node v = *(candidates->begin());
		nodesState->find(v)->second = label;
		community->insert(v);
		isCore = this->isCore(v);
		if (isCore.first) {
			for (node x : isCore.second) {
				int tmp = nodesState->find(x)->second;
				if (tmp < 0) {
					nodesState->find(x)->second = label;
					community->insert(x);
				}
				if (tmp == -1) {
					candidates->insert(x);
				}
			}
		}
		candidates->erase(candidates->find(v));
	}
}

} /* namespace NetworKit */
