/*
 * Partition.cpp
 *
 *  Created on: 03.10.2013
 *      Author: cls
 */

#include "Partition.h"

namespace NetworKit {

Partition::Partition(index z) : z(z), omega(0), data(z, none) { // FIXME: data is initialized with z+1, is that correct?

}

void Partition::addToSubset(index s, index e) {
	assert (data[e] == none);	// guarantee that element was unasigned
	assert (s <= omega); // do not create new subset ids
	data[e] = s;
}

void Partition::moveToSubset(index s, index e) {
	assert (s <= omega); // do not create new subset ids
	data[e] = s;
}

void Partition::toSingleton(index e) {
	data[e] = newSubsetId();
}

void Partition::allToSingletons() {
	for (index e = 0; e <= this->z; ++e) {
		toSingleton(e);
	}
}

void Partition::mergeSubsets(index s, index t) {
	assert (s <= omega);
	assert (t <= omega);
	index m = newSubsetId(); // new id for merged set
	for (index e = 0; e <= this->z; ++e) {
		if (data[e] == s || data[e] == t) {
			data[e] = m;
		}
	}
}

bool Partition::isOnePartition(const std::set<index>& elements) { //FIXME what for is elements needed?
	index one = data[0];	// first subset id should be equal to all others
	// TODO: use iterator forEntries and pair-wise comparison?
	for (index e = 0; e < this->z; ++e) { // FIXME constructor initializes data with z+1, so <= is necessary. 
		if (data[e] != one) {
			return false;
		}
	}
	return true;
}

bool Partition::isSingletonPartition(const std::set<index>& elements) const { //FIXME what for is elements needed?
	return (numberOfElements() == numberOfSubsets());
}

count Partition::numberOfSubsets() const {
	std::vector<int> exists(upperBound(), 0); // a boolean vector would not be thread-safe

	this->parallelForEntries([&](index e, index s) {
		if (s != none) {
			exists[s] = 1;
		}
	});

	count k = 0; // number of actually existing clusters
	#pragma omp parallel for reduction(+:k)
	for (index i = 0; i < upperBound(); ++i) {
		if (exists[i]) {
			k++;
		}
	}

	return k;
}

index Partition::upperBound() const {
	return omega + 1;
}

index Partition::lowerBound() const {
	return 0;
}

void Partition::compact() {
	// TODO: compact partition
	// PROBLEM: conflict between old and new IDs?
	std::vector<index> assignedSubsetIDs;
	std::map<index,index> mapToCompactedIDs; // create map to
	index i = 1;
	this->forEntries([&](index e, index s){ // get assigned SubsetIDs and create a map with new IDs
		if (mapToCompactedIDs.find(s) == mapToCompactedIDs.end()) {
			//std::find(assignedSubsetIDs.begin(),assignedSubsetIDs.end()) != assignedSubsetIDs.end()
			//((assignedSubsetIDs.push_back(s);
			mapToCompactedIDs[s] = i++;
		}
	});
/*	std::map<index,index> mapToCompactedIDs; // create map to
	for (index i = 0; i < assignedSubsetIDs.size();i++) {
		mapToCompactedIDs[assignedSubsetIDs[i]] = i;
	}*/
	this->forEntries([&](index e, index s){ // replace old SubsetIDs with the new IDs
		data[e] = mapToCompactedIDs[s];
	});
	omega = i-1;
	// TODO: test implementation
	
}

bool Partition::contains(index e) const {
	return (e <= z) && (data[e] != none);	// e is in the element index range and the entry is not empty
}

bool Partition::inSameSubset(index e1, index e2) const {
	assert (data[e1] != none);
	assert (data[e2] != none);
	return (data[e1] == data[e2]);
}

std::vector<count> Partition::subsetSizes() const {
	std::vector<count> sizes;
	std::map<index, count> map = this->subsetSizeMap();
	for (auto kv : map) {
		sizes.push_back(kv.second);
	}
	return sizes;
}

std::map<index, count> Partition::subsetSizeMap() const {
	std::map<index, count> subset2size;

	this->forEntries([&](index e, index s){
		if (s != none) {
			subset2size[s] += 1;
		}
	});

	return subset2size;
}

std::set<index> Partition::getMembers(const index s) const {
	assert (s <= omega);
	std::set<index> subset;
	for (index e = 0; e <= this->z; ++e) {
		if (data[e] == s) {
			subset.insert(e);
		}
	}
	return subset;
}

count Partition::numberOfElements() const {
	return z;	// z is the maximum element id
}


} /* namespace NetworKit */
