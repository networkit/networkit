/*
 * Partition.cpp
 *
 *  Created on: 03.10.2013
 *      Author: cls
 */

#include "Partition.h"

namespace NetworKit {

Partition::Partition() : z(0), omega(0), data(0) { 

}


Partition::Partition(index z) : z(z), omega(0), data(z, none) {  //z(z-1);data(z,none);

}

Partition::Partition(index z, index defaultValue) : z(z), omega(0), data(z, defaultValue) {  //z(z-1);data(z,none);

}

void Partition::addToSubset(index s, index e) {
	assert (data[e] == none);	// guarantee that element was unassigned
	assert (s <= omega);		// do not create new subset ids
	data[e] = s;
}

void Partition::moveToSubset(index s, index e) {
	assert (this->contains(e));
	assert (s <= omega); 		// do not create new subset ids
	data[e] = s;
}

void Partition::toSingleton(index e) {
	data[e] = newSubsetId();
}

void Partition::allToSingletons() {
	for (index e = 0; e < this->z; ++e) {
		toSingleton(e);
	}
}

index Partition::mergeSubsets(index s, index t) {
	assert (s <= omega);
	assert (t <= omega);
	if (s != t) {
		index m = newSubsetId(); // new id for merged set
		for (index e = 0; e < this->z; ++e) {
			if (data[e] == s || data[e] == t) {
				data[e] = m;
			}
		}
		return m;
	}
	return none; // no new cluster formed
}
/*
bool Partition::isOnePartition(Graph& G) { //FIXME what for is elements needed? const std::set<index>& elements
	index one = data[0];	// first subset id should be equal to all others
	// TODO: use iterator forEntries and pair-wise comparison?
	for (index e = 0; e < this->z; ++e) { // FIXME constructor initializes data with z+1, so <= is necessary. 
		if (data[e] != one) {
			return false;
		}
	}
	return true;
}*/

/*bool Partition::isSingletonPartition(Graph& G) const { //FIXME what for is elements needed? const std::set<index>& elements
	return (numberOfElements() == numberOfSubsets());
}
*/
index Partition::extend() {
	data.push_back(none);
	z++;
	assert (z == data.size()); //(data.size() - 1)
	return z-1;
}

void Partition::remove(index e) {
	assert (e < z);
	data[e] = none;
}

count Partition::numberOfSubsets() const {
	auto n = upperBound();
	std::vector<int> exists(n, 0); // a boolean vector would not be thread-safe
	this->parallelForEntries([&](index e, index s) {
		if (s != none) {
			exists[s] = 1;
		}
	});
	count k = 0; // number of actually existing clusters
	#pragma omp parallel for reduction(+:k)
	for (index i = 0; i < n; ++i) {
		if (exists[i]) {
			k++;
		}
	}
	return k;
}

void Partition::setUpperBound(index upper) {
	this->omega = upper-1;
}

index Partition::upperBound() const {
	return omega+1;
}

index Partition::lowerBound() const {
	return 0;
}

void Partition::compact() {
	std::map<index,index> compactingMap; // first index is the old partition index, "value" is the index of the compacted index
	index i = 0;
	this->forEntries([&](index e, index s){ // get assigned SubsetIDs and create a map with new IDs
		if (s!= none) { 
			auto result = compactingMap.insert(std::make_pair(s,i));
			if (result.second) ++i;
		}
	});
	this->parallelForEntries([&](index e, index s){ // replace old SubsetIDs with the new IDs
		data[e] = compactingMap[s];
	});
	this->setUpperBound(i); // does i contain the right value?
}

bool Partition::contains(index e) const {
	return (e < z) && (data[e] != none);	// e is in the element index range and the entry is not empty
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
	for (index e = 0; e < this->z; ++e) {
		if (data[e] == s) {
			subset.insert(e);
		}
	}
	return subset;
}

count Partition::numberOfElements() const {
	return z;	// z is the maximum element id
}

std::vector<index> Partition::getVector() {
	return this->data; //FIXME is this appropriate? - why not?
}


std::set<std::set<index> > Partition::getSubsets() {
	std::vector<std::set<index> > table(omega);
	this->forEntries([&](index e, index s){
		table[s].insert(e);
	});

	std::set<std::set<index> > subsets;
	for (auto set : table) {
		if (set.size() > 0) {
			subsets.insert(set);
		}
	}
	return subsets;
}

void Partition::allToOnePartition() {
	omega = 0;
	this->parallelForEntries([&](index e, index s) {
		this->data[e] = 0;
	});
}

std::string Partition::getName() const {
	return this->name;
}

void Partition::setName(std::string name) {
	this->name = name;
}

} /* namespace NetworKit */
