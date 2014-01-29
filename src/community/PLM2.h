/*
 * MLPLM.h
 *
 *  Created on: 20.11.2013
 *      Author: cls
 */

#ifndef PLMR2_H_
#define PLMR2_H_

#include "Clusterer.h"

namespace NetworKit {

/**
 * MultiLevel Parallel LocalMover - a multi-level modularity maximizer.
 */
class PLM2: public NetworKit::Clusterer {

public:

	/**
	 * @param[in]	refine	add a second move phase to refine the communities
	 * @param[in]	par		parallelization strategy
	 * @param[in]	gamma	multi-resolution modularity parameter:
	 * 							1.0 -> standard modularity
	 * 							0.0 -> one community
	 * 							2m 	-> singleton communities
	 * @param[in]	maxIter		maximum number of iterations for move phase	
	 *
	 */
	PLM2(bool refine=false, double gamma = 1.0, std::string par="balanced", count maxIter=50);


	std::string toString() const override;

	/**
	 * Detect communities in the given graph.
	 */
	Clustering run(Graph& G) override;

	static std::pair<Graph, std::vector<node>> coarsen(const Graph& G, const Clustering& zeta);

	static Clustering prolong(const Graph& Gcoarse, const Clustering& zetaCoarse, const Graph& Gfine, std::vector<node> nodeToMetaNode);

private:

	std::string parallelism;
	bool refine;
	double gamma = 1.0;
	count maxIter;
};

} /* namespace NetworKit */

#endif /* MLPLM_H_ */
