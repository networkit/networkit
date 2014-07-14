/*
 * Louvain.h
 *
 *  Created on: 25.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#ifndef PLMOLD_H_
#define PLMOLD_H_

#include "CommunityDetectionAlgorithm.h"

namespace NetworKit {

/**
 * @ingroup community
 * PLM - Parallel Louvain Method community detection algorithm
 * The Lovain method is a locally greedy procedure for maximizing modularity.
 * This is a parallel implementation.
 */
class PLMOld: public NetworKit::CommunityDetectionAlgorithm {


public:

	std::string VERSION;	// algorithm version number - increment in constructor for significant changes to the implementation

	/**
	 * @param[in]	par		parallelization strategy
	 * @param[in]	gamma	multi-resolution modularity parameter:
	 * 							1.0 -> standard modularity
	 * 							0.0 -> one community
	 * 							2m 	-> singleton communities
	 *
	 */
	PLMOld(std::string par="balanced", double gamma = 1.0);

	virtual Partition pass(const Graph& G);

	virtual Partition run(const Graph& G);

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;


protected:

	bool anyChange;	//!< indicates whether any change was made to the clustering in the last pass over the nodes
	std::string parallelism; //!< switch for the kind of parallelization strategy to use
	double gamma;	//!< multi-resolution modularity parameter
};

} /* namespace NetworKit */
#endif /* LOUVAIN_H_ */
