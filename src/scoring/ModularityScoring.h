/*
 * Modularity.h
 *
 *  Created on: 15.10.2012
 *      Author: cls
 */

#ifndef MODULARITY_H_
#define MODULARITY_H_

#include "EdgeScoring.h"

#include "../clustering/Clustering.h"


namespace EnsembleClustering {

// TODO: implement modularity as in Python prototype

class ModularityScoring: public EnsembleClustering::EdgeScoring {

public:

	ModularityScoring();

	virtual ~ModularityScoring();


	/**
	 * Returns an edge score for an edge (u,v) which expresses the
	 * modularity increase which can be gained by merging
	 * the clusters of u and v.
	 *
	 * @param[in]	u	source node id
	 * @param[out]	v	target node id
	 *
	 */
	virtual double scoreEdge(node u, node v);

	/**
	 * Calculates the modularity of the given clustering;
	 *
	 */
	virtual double mod(Clustering& clustering) =0; // TODO: needed?

	/**
	 * Calculates the difference in modularity that would result from a merger of
	 * two clusters.
	 *
	 */
	virtual double deltaMod(cluster c, cluster d) =0;

	virtual double cutweight(cluster c, cluster d) =0;

	virtual double weight(cluster c) =0;
};

} /* namespace EnsembleClustering */
#endif /* MODULARITY_H_ */
