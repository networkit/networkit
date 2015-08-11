/*
 * BloomFilter.cpp
 *
 *  Created on: 08.08.2015
 *      Author: Henning
 */

#include "BloomFilter.h"
#include "Random.h"

namespace Aux {

BloomFilter::BloomFilter(count numHashes, count size): numHashes(numHashes), size(size) {
	membership.resize(numHashes);
	salts.resize(numHashes);

	for (index i = 0; i < numHashes; ++i) {
		membership[i].resize(size, false);
		salts[i] = Aux::Random::integer();
	}
}

void BloomFilter::insert(index key) {
	for (index func = 0; func < numHashes; ++func) {
		membership[func][hash(key, func)] = true;
	}
}

bool BloomFilter::isMember(index key) const {
	bool result = true;
	for (index func = 0; func < numHashes; ++func) {
		result = result && membership[func][hash(key, func)];
	}
	return result;
}

index BloomFilter::hash(index key, index hfunc) const {
	index result = 0;
	std::hash<index> myhash;

	result = myhash(key ^ salts[hfunc]);
	return (result % size);
}

}
