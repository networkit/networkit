/*
 * SelectiveSCAN.cpp
 *
 *  Created on: 14.06.2013
 *      Author: cls
 */

#include "SelectiveSCAN.h"

namespace NetworKit {

SelectiveSCAN::SelectiveSCAN(): epsilon(0.7), nhu(3) {

}

SelectiveSCAN::~SelectiveSCAN() {
	// TODO Auto-generated destructor stub
}

std::unordered_map<node, std::unordered_set<node>> SelectiveSCAN::seedSetExpansion(
			Graph& G, std::vector<node> set){

	std::unordered_map<node, node> nodesState;
	std::unordered_map<node, std::unordered_set<node>> communitites;
	SelectiveSCAN SSCAN;

	for (node u : set) {
		nodesState.insert(std::pair<node,int>(u, -1));
	}
	for (node u : set) {
		if ((nodesState.find(u))->second == -1  && SSCAN.isCore(u, G).first){

		}
	}
	return communitites;
}

double SelectiveSCAN::nodeNodeSimilarity (node u, node v, Graph& G) {

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
	return (double) (inter + 2)/ (double) (uni +2);
}

std::pair<bool,std::vector<node>> SelectiveSCAN::isCore(node u, Graph& G) {

	bool core = false;
	SelectiveSCAN SSCAN;
	std::vector<node> similarNeighbors;
	int count = 0;
	G.forNeighborsOf(u, [&](node v){
		if (SSCAN.nodeNodeSimilarity(u, v, G) >= SSCAN.epsilon) {
			count++;
			similarNeighbors.push_back(v);
		}
	});
	if (count >= SSCAN.nhu) {
		core = true;
	}
	return std::pair<bool,std::vector<node>>(core, similarNeighbors);
}

} /* namespace NetworKit */
