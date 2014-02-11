/*
 * ClusterContracter.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERCONTRACTER_H_
#define CLUSTERCONTRACTER_H_


#include "GraphCoarsening.h"
#include "../clustering/Clustering.h"

namespace NetworKit {

class ClusterContracter: public GraphCoarsening {

public:

	ClusterContracter();

	virtual ~ClusterContracter();

	virtual std::pair<Graph, std::vector<node> > run(const Graph& G, const Clustering& zeta);




};


} // namespace

#endif /* CLUSTERCONTRACTER_H_ */
