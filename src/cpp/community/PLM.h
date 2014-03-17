/*
 * MLPLM.h
 *
 *  Created on: 20.11.2013
 *      Author: cls
 */

#ifndef PLM_H_
#define PLM_H_

#include "CommunityDetectionAlgorithm.h"

namespace NetworKit {

/**
 * MultiLevel Parallel LocalMover - a multi-level modularity maximizer.
 */
class PLM: public NetworKit::CommunityDetectionAlgorithm {

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
	PLM(bool refine=false, double gamma = 1.0, std::string par="balanced", count maxIter=32);


	std::string toString() const override;

	/**
	 * Detect communities in the given graph.
	 */
	Partition run(Graph& G) override;

	static std::pair<Graph, std::vector<node>> coarsen(const Graph& G, const Partition& zeta);

	static Partition prolong(const Graph& Gcoarse, const Partition& zetaCoarse, const Graph& Gfine, std::vector<node> nodeToMetaNode);

private:

	std::string parallelism;
	bool refine;
	double gamma = 1.0;
	count maxIter;
};

} /* namespace NetworKit */

#endif /* MLPLM_H_ */
