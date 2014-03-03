/*
 * Louvain.h
 *
 *  Created on: 25.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#ifndef LOUVAIN_H_
#define LOUVAIN_H_

#include "Clusterer.h"

namespace NetworKit {

/**
 * PLM - Parallel Louvain Method community detection algorithm
 * The Lovain method is a locally greedy procedure for maximizing modularity. 
 * This is a parallel implementation.
 */
class PLM: public NetworKit::Clusterer {


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
	PLM(std::string par="balanced", double gamma = 1.0);

	virtual ~PLM();

	virtual Partition pass(Graph& G);

	virtual Partition run(Graph& G);

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
