/*
 * UnionFind.cpp
 *
 *  Created on: 11.08.2014
 *      Author: Marcel Radermacher
 */

#include "UnionFind.h"

namespace NetworKit {


void UnionFind::allToSingletons() {
	for (index i = 0; i < data.size(); ++i) {
		data[i] = -1;
	}
}

index UnionFind::find(index u) {
	if (data[u] >= 0) {
		data[u] = find(data[u]);
		return data[u];
	} else {
		return u;
	}
}

void UnionFind::merge(index u, index v) {
	index set_u = find(u);
	index set_v = find(v);
	if (set_u == set_v) return;

	if (data[set_u] > data[set_v]) {
		data[set_u] = set_v;
	} else if (data[set_v] > data[set_u]) {
		data[set_v] = set_u;
	} else {
		data[set_u] = set_v;
		--data[set_v];
	}
}

Partition UnionFind::toPartition() {
	Partition p(data.size());
	p.setUpperBound(data.size());
	for (index e = 0; e < data.size(); ++e) {
		p.moveToSubset(find(e), e);
	}	
	return std::move(p);
}


}
