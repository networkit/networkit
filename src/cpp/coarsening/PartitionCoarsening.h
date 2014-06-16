/*
 * PartitionCoarsening.h
 *
 *  Created on: 28.01.2014
 *      Author: cls
 */

#ifndef PARTITIONCOARSENING_H_
#define PARTITIONCOARSENING_H_

#include "GraphCoarsening.h"
#include "../structures/Partition.h"

namespace NetworKit {

/*
 *
 */
class PartitionCoarsening: public NetworKit::GraphCoarsening {

public:

	virtual std::pair<Graph, std::vector<node> > run(const Graph& G, const Partition& zeta);



};

} /* namespace NetworKit */

#endif /* PARTITIONCOARSENING_H_ */
