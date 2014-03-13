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
		int neighborhood1 = G.degree(u);
		int neighborhood2 = G.degree(v);

		G.forNeighborsOf(u, [&](node x){
			if (x != v && x!= u ) {
				if (G.hasEdge(x, v)) {
					inter++;
				}
			}
		});
		if (!G.hasEdge(u, u)) {
			neighborhood1++;
		}
		if (!G.hasEdge(v, v)) {
			neighborhood2++;
		}

		return (1 - ((double) (inter + 2)) / (sqrt(neighborhood2 * neighborhood1)));
}

} /* namespace NetworKit */
