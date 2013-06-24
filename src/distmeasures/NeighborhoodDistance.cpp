/*
 * NeighborhoodDistance.cpp
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#include "NeighborhoodDistance.h"

namespace NetworKit {

NeighborhoodDistance::NeighborhoodDistance(const Graph& G) : NodeDistance(G) {
	// TODO Auto-generated constructor stub

}

NeighborhoodDistance::~NeighborhoodDistance() {
	// TODO Auto-generated destructor stub
}

void NeighborhoodDistance::preprocess() {
	// no preprocessing necessary
}

double NeighborhoodDistance::distance(node u, node v) {
	count inter = 0;
	count uni = 0;
	G.forNeighborsOf(u, [&](node x){
		if (x != v && x!= u ) {
			if (G.hasEdge(x, v)) {
				inter++;
				uni++;
			} else {
				uni++;
			}
		}
	});
	G.forNeighborsOf(v, [&](node x){
		if (x != u && x != v) {
			if (!G.hasEdge(x, u)) {
				uni++;
			}
		}
	});
	return 1 - ((double) (inter + 2) / (double) (uni +2));
}

} /* namespace NetworKit */
