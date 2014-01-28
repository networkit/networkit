/*
 * PartitionCoarsening.h
 *
 *  Created on: 28.01.2014
 *      Author: cls
 */

#ifndef PARTITIONCOARSENING_H_
#define PARTITIONCOARSENING_H_

#include "GraphCoarsening.h"
#include "../clustering/Clustering.h"

namespace NetworKit {

/*
 *
 */
class PartitionCoarsening: public NetworKit::GraphCoarsening {

public:

	virtual std::pair<Graph, std::vector<node> > run(Graph& G, Clustering& zeta);



};

} /* namespace NetworKit */

#endif /* PARTITIONCOARSENING_H_ */
