/*
 * HashingOverlapper.cpp
 *
 *  Created on: 31.01.2013
 *      Author: cls
 */

#include "HashingOverlapper.h"

namespace EnsembleClustering {

HashingOverlapper::HashingOverlapper() {
	// TODO Auto-generated constructor stub

}

HashingOverlapper::~HashingOverlapper() {
	// TODO Auto-generated destructor stub
}

Clustering HashingOverlapper::run(Graph& G,
		std::vector<Clustering>& clusterings) {

	// hash function sdbm
	auto sdbm = [](int64_t cid) {
		unsigned char* str = (unsigned char*) &cid;
		unsigned long h = 0;
		int c;
		while (c = *str++) {
			h = c + (h << 6) + (h << 16) - h;
		}
		return h;
	};

	auto djb2 = [](int64_t cid) {
		unsigned char* str = (unsigned char*) &cid;
		unsigned long hash = 5381;
		int c;
		while (c = *str++) {
			hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
		}
		return hash;
	};

	Clustering core(G.numberOfNodes());

	// select hash function
	auto hash = djb2;

	core.setAll(0);
	for (Clustering zeta : clusterings) {
		G.parallelForNodes([&](node v){
			core[v] += (hash(zeta[v]) & 0xffff);
		});
	}

	core.compact();


//	// TODO: this is a quick fix for the issue that clusterings start with an upper id bound of n. bounds are needed and checked for iteration over all clusters
//	cluster maxCluster = 0;
//	G.forNodes([&](node v) {
//		if (core[v] > maxCluster) {
//			maxCluster = core[v];
//		}
//	});
//	INFO("upperBound: " << core.upperBound());
//	core.setUpperBound(maxCluster + 1);
//	INFO("upperBound: " << core.upperBound());
//	assert (core.upperBound() == (maxCluster + 1));
//	INFO("maxCluster: " << maxCluster);

	return core;
}

} /* namespace EnsembleClustering */
