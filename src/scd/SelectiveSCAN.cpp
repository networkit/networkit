/*
 * SelectiveSCAN.cpp
 *
 *  Created on: 14.06.2013
 *      Author: cls
 */

#include "SelectiveSCAN.h"

namespace NetworKit {

SelectiveSCAN::SelectiveSCAN(Graph& G): SelectiveCommunityDetector(G), epsilon(0.7), mu(3) {

}

SelectiveSCAN::~SelectiveSCAN() {
	// TODO Auto-generated destructor stub
}

std::unordered_map<node, std::unordered_set<node>> SelectiveSCAN::run(std::unordered_set<node> set){

	std::unordered_map<node, node> nodesState;
	std::unordered_map<node, std::unordered_set<node>> communitites;

	for (node u : set) {
		nodesState.insert(std::pair<node,int>(u, -1));
	}
	for (node u : set) {
		if ((nodesState.find(u))->second == -1  && this->isCore(u).first){

		}
	}
	return communitites;
}

double SelectiveSCAN::nodeDistance(node u, node v) {

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

std::pair<bool,std::vector<node>> SelectiveSCAN::isCore(node u) {

	bool core = false;
	std::vector<node> similarNeighbors;
	int count = 0;
	G.forNeighborsOf(u, [&](node v){
		if (this->nodeDistance(u, v) <= this->epsilon) {
			count++;
			similarNeighbors.push_back(v);
		}
	});
	if (count >= this->mu) {
		core = true;
	}
	return std::pair<bool,std::vector<node>>(core, similarNeighbors);
}

} /* namespace NetworKit */
