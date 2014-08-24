/*
 * Overlapper.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef OVERLAPPER_H_
#define OVERLAPPER_H_

#include <set>
#include <vector>

#include "../graph/Graph.h"
#include "../structures/Partition.h"

namespace NetworKit {


/**
 * @ingroup overlap
 * Abstract base class for algorithms which determine the overlap of multiple partitions.
 */
class Overlapper {

public:

	virtual Partition run(const Graph& G, const  std::vector<Partition>& clusterings) = 0;

};

} /* namespace NetworKit */
#endif /* OVERLAPPER_H_ */
