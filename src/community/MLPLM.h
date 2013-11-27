/*
 * MLPLM.h
 *
 *  Created on: 20.11.2013
 *      Author: cls
 */

#ifndef MLPLM_H_
#define MLPLM_H_

#include "Clusterer.h"

namespace NetworKit {

/**
 * MultiLevel Parallel LocalMover - a multi-level modularity maximizer.
 */
class MLPLM: public NetworKit::Clusterer {

public:

	/**
	 * @param[in]	par		parallelization strategy
	 * @param[in]	refine	add a second move phase to refine the communities
	 * @param[in]	gamma	multi-resolution modularity parameter:
	 * 							1.0 -> standard modularity
	 * 							0.0 -> one community
	 * 							2m 	-> singleton communities
	 *
	 */
	MLPLM(std::string par="simple",  bool refine=true, double gamma = 1.0);


	std::string toString() const override;

	/**
	 * Detect communities in the given graph.
	 */
	Clustering run(Graph& G) override;

	std::pair<Graph, std::vector<node>> coarsen(const Graph& G, const Clustering& zeta);

	Clustering prolong(const Graph& Gcoarse, const Clustering& zetaCoarse, const Graph& Gfine, std::vector<node> nodeToMetaNode);

private:

	std::string parallelism;
	bool refine;
	double gamma = 1.0;
};

} /* namespace NetworKit */

#endif /* MLPLM_H_ */
