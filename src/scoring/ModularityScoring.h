/*
 * Modularity.h
 *
 *  Created on: 15.10.2012
 *      Author: cls
 */

#ifndef MODULARITY_H_
#define MODULARITY_H_

#include "EdgeScoring.h"

#include "../clustering/base/Clustering.h"


namespace EnsembleClustering {

// TODO: implement modularity as in Python prototype

class ModularityScoring: public EnsembleClustering::EdgeScoring {

protected:

	double omegaE;	//!< total weight of the graph

public:

	/**
	 * @param[in]	G	a graph instance
	 *
	 * Do not modify the graph while using this instance of ModularityScoring.
	 */
	ModularityScoring(Graph& G);

	virtual ~ModularityScoring();


	/**
	 * Returns an edge score for an edge (u,v) which expresses the
	 * modularity increase which can be gained by merging
	 * the clusters of u and v.
	 *
	 *		 $$\Delta mod(c, d) := \frac{1}{2 \omega(E)} \left ( 2 \omega(E) \omega(c,d) - \omega(c) \omega(d) \right ) $$
	 *
	 * @param[in]	u	source node id
	 * @param[out]	v	target node id
	 *
	 */
	virtual double scoreEdge(node u, node v);


//	/**
//	 * Calculates the difference in modularity that would result from a merger of
//	 * two clusters.
//	 *
//	 */
//	virtual double deltaMod(cluster c, cluster d) =0;
//
//	virtual double cutweight(cluster c, cluster d) =0;
//
//	virtual double weight(cluster c) =0;
};

} /* namespace EnsembleClustering */
#endif /* MODULARITY_H_ */
