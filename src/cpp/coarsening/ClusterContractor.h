/*
 * ClusterContractor.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERCONTRACTOR_H_
#define CLUSTERCONTRACTOR_H_


#include "GraphCoarsening.h"
#include "../structures/Partition.h"

namespace NetworKit {

class ClusterContractor: public GraphCoarsening {

public:

	ClusterContractor();

	virtual ~ClusterContractor();

	virtual std::pair<Graph, std::vector<node> > run(const Graph& G, const Partition& zeta);




};


} // namespace

#endif /* CLUSTERCONTRACTOR_H_ */
