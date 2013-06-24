/*
 * TNeighborhoodDistance.cpp
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#include "TNeighborhoodDistance.h"

namespace NetworKit {

TNeighborhoodDistance::TNeighborhoodDistance(const Graph& G) : TNodeDistance(G) {
	// TODO Auto-generated constructor stub

}

TNeighborhoodDistance::~TNeighborhoodDistance() {
	// TODO Auto-generated destructor stub
}

void TNeighborhoodDistance::initialize(const Parameters& param) {
	// no initialization needed
}

double TNeighborhoodDistance::distance(node u, node v) {
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
