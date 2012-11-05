/*
 * Modularity.h
 *
 *  Created on: 15.10.2012
 *      Author: cls
 */

#ifndef MODULARITY_H_
#define MODULARITY_H_

#include "EdgeScoring.h"
#include "../graph/Graph.h"


namespace EnsembleClustering {


// TODO: implement modularity as in Python prototype

class Modularity: public EnsembleClustering::EdgeScoring {

public:

	Modularity();

	virtual ~Modularity();


	/**
	 * Returns an edge score for an edge (u,v) which expresses the
	 * modularity increase which can be gained by merging
	 * the clusters of u and v.
	 *
	 * @param[in]	u	source node id
	 * @param[out]	v	target node id
	 *
	 */
	virtual double scoreEdge(Node u, Node v);

	/**
	 * Calculates the modularity of the given clustering;
	 *
	 */
	double mod(Clustering* clustering);

	/**
	 * Calculates the difference in modularity that would result from a merger of
	 * two clusters.
	 *
	 */
	virtual double deltaMod(Cluster* c, Cluster* d);

	double cutweight(Cluster* c, Cluster* d);

	double weight(Cluster* c);
};

} /* namespace EnsembleClustering */
#endif /* MODULARITY_H_ */
