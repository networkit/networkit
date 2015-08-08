/*
 * BloomFilter.cpp
 *
 *  Created on: 08.08.2015
 *      Author: Henning
 */

#include "BloomFilter.h"
#include "Random.h"

namespace Aux {

BloomFilter::BloomFilter(count numHashes): numHashes(numHashes) {
	size = 6291469;
	membership.resize(numHashes);
	salts.resize(numHashes);

	for (index i = 0; i < numHashes; ++i) {
		membership[i].resize(size, false);
		salts[i] = Aux::Random::integer();
	}

//	DEBUG("end of constructor");
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
//	TRACE("result: ", result);
	return (result % size);
}

}
