/*
 * Cover.cpp
 *
 *  Created on: 03.10.2013
 *      Author: cls
 */

#include "Cover.h"

#include <algorithm>

namespace NetworKit {

Cover::Cover(index z) : z(z), omega(0), data(z+1) {
}

bool Cover::contains(index e) const {
	return (e <= z) && (! data[e].empty());	// e is in the element index range and the entry is not empty
}

bool Cover::inSameSubset(index e1, index e2) const {
	assert (e1 <= z);
	assert (e2 <= z);
	assert (! data[e1].empty());
	assert (! data[e2].empty()); // elements cannot be unassigned - it may be possible to change this behavior
	std::set<index> intersect;
	std::set_intersection(data[e1].begin(),data[e1].end(),data[e2].begin(), data[e2].end(), std::inserter(intersect,intersect.begin()));
	return (!intersect.empty());
}

std::set<index> Cover::getMembers(const index s) const {
	assert (s <= omega);
	std::set<index> members;
	for (index e = 0; e <= this->z; ++e) {
		for (index t : data[e]) {
			if (t == s) {
				members.insert(e);
			}
		}
	}
	return members;
}

void Cover::addToSubset(index s, index e) {
	assert (e <= z);
	assert (s <= omega);
	data[e].insert(s);
}

void Cover::moveToSubset(index s, index e) {
	assert (e <= z);
	assert (s <= omega);
	data[e].clear();
	data[e].insert(s);
}

void Cover::toSingleton(index e) {
	assert (e <= z);
	data[e].clear();
	data[e].insert(newSubsetId());
}

void Cover::allToSingletons() {
	// TODO:
}

void Cover::mergeSubets(index s, index t) {
	// TODO:
}

index Cover::upperBound() const {
	// TODO:
}

index Cover::lowerBound() const {
	// TODO:
}

std::vector<count> Cover::subsetSizes() const {
	// TODO:
}

std::map<index, count> Cover::subsetSizeMap() const {
	// TODO:
}


} /* namespace NetworKit */

