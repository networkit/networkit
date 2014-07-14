/*
 * ParallelPartitionCoarsening.h
 *
 *  Created on: 03.07.2014
 *      Author: cls
 */

#ifndef PARALLELPARTITIONCOARSENING_H_
#define PARALLELPARTITIONCOARSENING_H_

#include "../Globals.h"
#include "GraphCoarsening.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * @ingroup coarsening
 */
class ParallelPartitionCoarsening: public NetworKit::GraphCoarsening {

public:

	virtual std::pair<Graph, std::vector<node> > run(const Graph& G, const Partition& zeta);


private:

	class ThreadSafeAdjacency {

	public:

		ThreadSafeAdjacency(count z, count nThreads) : nThreads(nThreads), z(z), adja(z), weight(z) {
			for (count i = 0; i < z; ++i) {
				adja[i].resize(nThreads);
				weight[i].resize(nThreads);
			}
		};

	private:
		count nThreads;
		count z;
		std::vector< std::vector<std::vector<node> > > adja;
		std::vector< std::vector<std::vector<edgeweight> > > weight;
	};

};

} /* namespace NetworKit */

#endif /* PARALLELPARTITIONCOARSENING_H_ */
