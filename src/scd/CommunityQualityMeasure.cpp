/*
 * CommunityQualityMeasure.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "CommunityQualityMeasure.h"

namespace NetworKit {

CommunityQualityMeasure::CommunityQualityMeasure() {
	// TODO Auto-generated constructor stub

}

CommunityQualityMeasure::~CommunityQualityMeasure() {
	// TODO Auto-generated destructor stub
}



LocalModularity::LocalModularity() {
}


LocalModularity::~LocalModularity() {
}

double LocalModularity::getQuality(
		const std::unordered_set<node>& community, const Graph& G) {
	double inside = 0;
	double outside = 0;
	std::unordered_set<node> boundary;

	for (node u : community) {
		G.forNeighborsOf(u, [&](node x){
			if (community.find(x) == community.end()){
				outside ++;
				if (boundary.find(u) == boundary.end()) {
					boundary.insert(u);
				}
			} else {
				if (u == x) {
					inside++;
				} else {
					inside = inside + 0.5;
				}
			}
		});
	}
	return (inside / community.size()) / (outside / boundary.size());
}


} /* namespace NetworKit */
