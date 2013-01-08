/*
 * ClusteringProjector.h
 *
 *  Created on: 07.01.2013
 *      Author: cls
 */

#ifndef CLUSTERINGPROJECTOR_H_
#define CLUSTERINGPROJECTOR_H_

#include "../clustering/Clustering.h"
#include "GraphContraction.h"

namespace EnsembleClustering {

class ClusteringProjector {

public:

	ClusteringProjector();

	virtual ~ClusteringProjector();

	/**
	 * Given
	 * 		@param[in]	contraction	a graph contraction, including a fine graph, a coarse graphm, and a mapping from fine to coarse
	 * 		@param[in]	zetaCoarse	a clustering of the coarse graph
	 * 	, project the clustering back to the fine graph to create a clustering of the fine graph.
	 * 		@param[out] 			a clustering of the fine graph
	 **/
	virtual Clustering projectBack(GraphContraction& contraction, Clustering& zetaCoarse);

};

} /* namespace EnsembleClustering */
#endif /* CLUSTERINGPROJECTOR_H_ */
