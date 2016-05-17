/*
 * HashingOverlapper.cpp
 *
 *  Created on: 31.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "HashingOverlapper.h"
#include "../auxiliary/Log.h"

namespace NetworKit {


Partition HashingOverlapper::run(const Graph& G, const std::vector<Partition>& clusterings) {

	DEBUG("Starting hashing overlapper");

	// hash function sdbm
	/*auto sdbm = [](int64_t cid) {
		unsigned char* str = (unsigned char*) &cid;
		unsigned long h = 0;
		int c;
		while ((c = *str++)) {
			h = c + (h << 6) + (h << 16) - h;
		}
		return h;
	};*/

	auto djb2 = [](int64_t cid) {
		unsigned char* str = (unsigned char*) &cid;
		unsigned long hash = 5381;
		int c;
		while ((c = *str++)) {
			hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
		}
		return hash;
	};

	Partition core(G.upperNodeIdBound());

	// select hash function
	auto hash = djb2;

	const count numC = clusterings.size();
	switch(numC) {
	case 0: {
		ERROR("No clustering provided! Will return 1-clustering!");
		core.allToOnePartition();
		return core;
	}
	case 1: {
		core = clusterings[0];
		break;
	}
	case 2: {
		const Partition& first = clusterings[0];
		const Partition& second = clusterings[1];

		// Assumption: second has at least as many nodes as first
		// TODO: accelerate if necessary and possible by moving if statements out of loop
		G.parallelForNodes([&](node v) {
			if (v >= first.numberOfElements()) {
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
		break;
	}
	default: {
		// Here we need to ensure that core is initialized to the 1-clustering with ID 0
		core.allToOnePartition();

		for (index c = 0; c < numC; ++c) {
			const Partition& zeta = clusterings[c];
			zeta.parallelForEntries([&](node v, index clv) {
				core[v] += (hash((c+2) * clv) & 0xffff);
			});
		}
	}
	}

	core.compact();

	return core;
}

} /* namespace NetworKit */
