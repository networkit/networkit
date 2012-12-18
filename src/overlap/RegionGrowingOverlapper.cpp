/*
 * RegionGrowingOverlapper.cpp
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#include "RegionGrowingOverlapper.h"

namespace EnsembleClustering {

RegionGrowingOverlapper::RegionGrowingOverlapper() {
	// TODO Auto-generated constructor stub

}

RegionGrowingOverlapper::~RegionGrowingOverlapper() {
	// TODO Auto-generated destructor stub
}

void RegionGrowingOverlapper::run(Graph& G, std::vector<Clustering> clusterings) {
	int64_t n = G.numberOfNodes();
	Clustering zetaOverlap(n);

	node r;
	G.breadthFirstEdgesFrom(r, [&](node u, node v) {
		bool together = true;
		for (Clustering& zeta : clusterings) {
			together = together && (zeta.clusterOf(u) == zeta.clusterOf(v));
		}
		if (together) {
			zetaOverlap.mergeClusters(zetaOverlap.clusterOf(u), zetaOverlap.clusterOf(v));
		}
	});

	// TODO: design and implement

}

} /* namespace EnsembleClustering */
