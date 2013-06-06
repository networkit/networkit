/*
 * CONductance.cpp
 *
 *  Created on: 06.06.2013
 *  AUthor Yassine Marrakchi
 */

#include "Conductance.h"


namespace NetworKit {


Conductance::Conductance() : LocalQualityMeasure() {
}

Conductance::~Conductance() {
	// TODO Auto-generated destructor stub
}

double Conductance::getQuality(std::unordered_set<node>& C, Graph& G) {

	double volume = 0;
	double boundary = 0;
	double all = 0;

	for (auto it = C.begin(); it != C.end(); ++it) {
		volume = volume + G.degree(*it);
		G.forNeighborsOf(*it, [&](node v){
			if (C.find(v) == C.end()) boundary++;
		});
	}

	G.forNodes([&](node v){
		all = all + G.degree(v);
	});

	if (volume == 0 || all-volume == 0)
		return 1;
	return boundary / std::min(volume, all-volume);
}

} /* namespace NetworKit */
