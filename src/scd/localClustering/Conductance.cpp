/*
 * CONductance.cpp
 *
 *  Created on: 06.06.2013
 *  AUthor Yassine Marrakchi
 */

#include "Conductance.h"
#include <algorithm>

namespace NetworKit {


Conductance::Conductance() : LocalQualityMeasure() {
}

Conductance::~Conductance() {
	// TODO Auto-generated destructor stub
}

double getQuality(std::unordered_set<node>& zeta, Graph& G) {

	double volume = 0;
	double boundary = 0;
	double all = 0;

	for (auto it = zeta.begin(); it != zeta.end(); ++it) {
		volume = volume + G.degree(*it);
		G.forNeighborsOf(*it, [&](node v){
			if (zeta.find(v) == zeta.end()) boundary++;
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
