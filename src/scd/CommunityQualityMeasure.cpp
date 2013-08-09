/*
 * CommunityQualityMeasure.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "CommunityQualityMeasure.h"

namespace NetworKit {

CommunityQualityMeasure::CommunityQualityMeasure(Graph& G) {
	this->G = &G;
	this->degSum = this->G->parallelSumForNodes([&](node u){
		return this->G->degree(u);
	});
}

CommunityQualityMeasure::~CommunityQualityMeasure() {
	// TODO Auto-generated destructor stub
}



LocalModularity::LocalModularity(Graph&G) : CommunityQualityMeasure(G) {
}


LocalModularity::~LocalModularity() {
}

double LocalModularity::getQuality(
	const std::unordered_set<node>& community) {
	double inside = 0;
	double outside = 0;
	std::unordered_set<node> boundary;

	for (node u : community) {
		this->G->forNeighborsOf(u, [&](node x){
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
	if (community.size() == 0) {
		return 0;
	}
	if(outside == 0 || boundary.size() == 0) {
		return inside / community.size();
	}
	return (inside / community.size()) / (outside / boundary.size());
}

Conduct::Conduct(Graph& G) : CommunityQualityMeasure(G) {
}

Conduct::~Conduct() {
}

double Conduct::getQuality(const std::unordered_set<node>& community) {
	count outside = 0;
	count volume = 0;
	for (node u : community) {
		this->G->forNeighborsOf(u, [&](node x){
			if (community.find(x) == community.end()){
				outside ++;
			}
		});
		volume = volume + G->degree(u);
	}
	if (community.empty()) {
		return 1;
	}
	if (volume == 0) {
		return 1;
	} else if (degSum-volume == 0) {
		return 1;
	}
	return ((double)outside)/ ((double)std::min((degSum-volume),volume));
}

} /* namespace NetworKit */


