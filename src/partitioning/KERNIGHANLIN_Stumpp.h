/*
 * Louvain.h
 *
 *  Created on: 25.02.2013
 *      Author: Matthias Stumpp
 */

#ifndef KERNIGHAN_LIN_H_
#define KERNIGHAN_LIN_H_

#include "../community/Clusterer.h"
#include "../clustering/Clustering.h"
#include "../clustering/ClusteringGenerator.h"
#include "../graph/Graph.h"
#include "../clustering/EdgeCut.h"

namespace NetworKit {

class KERNIGHAN_LIN: public NetworKit::Clusterer {

public:

	KERNIGHAN_LIN();

	virtual ~KERNIGHAN_LIN();

	virtual Clustering run(Graph& graph);

	virtual Clustering partition(Graph& g, Clustering& c);

	virtual void computeGainAllNodes(Graph& g, Clustering& c, NodeMap<double>& gain);

	virtual void computeGainForNode(Graph& g, node u, Clustering& c, NodeMap<double>& gain);

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;
};

} /* namespace NetworKit */
#endif /* KERNIGHAN_LIN_H_ */
