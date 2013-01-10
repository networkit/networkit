/*
 * ClusterContracter.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef CLUSTERCONTRACTER_H_
#define CLUSTERCONTRACTER_H_


#include "Contracter.h"
#include "../clustering/base/Clustering.h"
#include "../base/IndexMap.h"
#include "../coarsening/GraphContraction.h"

namespace EnsembleClustering {

class ClusterContracter: public Contracter {

public:

	ClusterContracter();

	virtual ~ClusterContracter();

	// TODO: return GraphContraction
	virtual GraphContraction run(Graph& G, Clustering& zeta);
};


} // namespace

#endif /* CLUSTERCONTRACTER_H_ */
