/*
 * BloomFilter.h
 *
 *  Created on: 08.08.2015
 *      Author: Henning
 */

#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_

#include "../auxiliary/Log.h"
#include "../Globals.h"

namespace Aux {

typedef NetworKit::index index;
typedef NetworKit::count count;

class BloomFilter {
private:
	count numHashes;
	count size;
	std::vector<std::vector<char> > membership;
	std::vector<index> salts;

	index hash(index key, index hfunc) const;

public:
	BloomFilter(count numHashes);
	virtual ~BloomFilter() = default;

	void insert(index key);

	bool isMember(index key) const;
};

}

#endif /* BLOOMFILTER_H_ */
