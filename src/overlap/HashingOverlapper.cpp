/*
 * HashingOverlapper.cpp
 *
 *  Created on: 31.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "HashingOverlapper.h"

namespace NetworKit {

HashingOverlapper::HashingOverlapper() {
	// TODO Auto-generated constructor stub

}

HashingOverlapper::~HashingOverlapper() {
	// TODO Auto-generated destructor stub
}

Clustering HashingOverlapper::run(Graph& G,
		std::vector<Clustering>& clusterings) {

	DEBUG("Starting hashing overlapper");

	// hash function sdbm
	auto sdbm = [](int64_t cid) {
		unsigned char* str = (unsigned char*) &cid;
		unsigned long h = 0;
		int c;
		while ((c = *str++)) {
			h = c + (h << 6) + (h << 16) - h;
		}
		return h;
	};

	auto djb2 = [](int64_t cid) {
		unsigned char* str = (unsigned char*) &cid;
		unsigned long hash = 5381;
		int c;
		while ((c = *str++)) {
			hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
		}
		return hash;
	};

	count n = G.numberOfNodes();
	Clustering core(n);

	// select hash function
	auto hash = djb2;

	core.setAll(0);
	const count numC = clusterings.size();
	const count upperId = clusterings.back().upperBound();
	const count summand = 341;
	if (numC > 2) {
		for (index c = 0; c < numC; ++c) {
			Clustering& zeta = clusterings[c];
			zeta.parallelForEntries([&](node v, cluster clv) {
				core[v] += (hash((c+2) * clv) & 0xffff);
			});
		}
	}
	else {
		Clustering& first = clusterings[0];
		Clustering& second = clusterings[1];

		// Assumption: second has at least as many nodes as first
		G.parallelForNodes([&](node v) {
			if (v >= first.numberOfEntries()) {
				core[v] = none;
			}
			else {
				if (first[v] == none || second[v] == none) {
					core[v] = none;
				}
				else {
					count key = ((first[v] ^ 0xffff) << 16) | (second[v] ^ 0xffff);
					core[v] = hash(key);
				}
			}
		});
	}

	core.compact();

	return core;
}

} /* namespace NetworKit */
