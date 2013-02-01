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


Clustering HashingOverlapper::run(Graph& G, std::vector<Clustering>& clusterings) {

	// hash function sdbm
	auto sdbm = [](int64_t cid) {
		unsigned char* str = (unsigned char*) &cid;
		unsigned long h = 0;
		int c;
		while (c = *str++)
		   h = c + (h << 6) + (h << 16) - h;
		return h;
	};


	Clustering core(G.numberOfNodes());

	// select hash function
	auto hash = sdbm;



	G.forallNodes([&](node v) {
		size_t cHash = 0;
		for (Clustering zeta : clusterings) {
			cHash += hash(zeta[v]);
		}
		cHash = cHash & 0xffff;
		core[v] = cHash;
	}, "parallel");

	// TODO: this is a quick fix for the issue that clusterings start with an upper id bound of n. bounds are needed and checked for iteration over all clusters
	cluster maxCluster = 0;
	G.forallNodes([&](node v){
		if (core[v] > maxCluster) {
			maxCluster = core[v];
		}
	});
	core.setUpperBound(maxCluster);

	return core;
}

} /* namespace EnsembleClustering */
