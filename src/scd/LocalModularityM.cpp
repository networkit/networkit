/*
 * LocalModularityM.cpp
 *
 *  Created on: 06.06.2013
 *  AUthor Yassine Marrakchi
 */

#include "LocalModularityM.h"


namespace NetworKit {


LocalModularityM::LocalModularityM() : LocalQualityMeasure() {
}

LocalModularityM::~LocalModularityM() {
	// TODO Auto-generated destructor stub
}

double LocalModularityM::getQuality(std::unordered_set<node>& zeta, Graph& G) {

	double inside = 0;
	double outside = 0;

	for (auto it = zeta.begin(); it != zeta.end(); ++it) {
		G.forNeighborsOf(*it, [&](node v){
			if (zeta.find(v) == zeta.end()){
				outside ++;
			} else {
				if (*it == v) {
					inside++;
				} else {
					inside = inside + 0.5;
				}
			}
		});
	}

	return inside / outside;
}

} /* namespace NetworKit */
