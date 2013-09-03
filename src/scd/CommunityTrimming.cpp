/*
 * CommunityTrimming.cpp
 *
 *  Created on: 10.06.2013
 *      Author: cls
 */

#include "CommunityTrimming.h"

namespace NetworKit {


NetworKit::CommunityTrimming::CommunityTrimming() {
}

NetworKit::CommunityTrimming::~CommunityTrimming() {
}

NetworKit::BoundarySharpness::BoundarySharpness():CommunityTrimming() {
}

NetworKit::BoundarySharpness::~BoundarySharpness() {
}

std::unordered_set<node> NetworKit::BoundarySharpness::run(
		std::unordered_set<node>& community, const Graph& G) {

	std::unordered_set<node> outliers;

	int in = 0;
	int out = 0;
	for (node u : community) {
		in = 0;
		out = 0;
		G.forNeighborsOf(u, [&](node v){
			if (community.find(v) == community.end()) {
				out++;
			} else {
				in++;
			}
		});
		if (out > in) {
			outliers.insert(u);
		}
	}

	for(node u : outliers) {
		community.erase(community.find(u));
	}
	return community;
}

NetworKit::DummyTrimming::DummyTrimming():CommunityTrimming() {
}

NetworKit::DummyTrimming::~DummyTrimming() {
}

std::unordered_set<node> NetworKit::DummyTrimming::run(
		std::unordered_set<node>& community, const Graph& G) {
	return community;
}

} /* namespace NetworKit */
