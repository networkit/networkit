/*
 * PLP.h
 *
 *  Created on: 07.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef PLP_H_
#define PLP_H_

#include "CommunityDetectionAlgorithm.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 *     As described in Ovelgoenne et al: An Ensemble Learning Strategy for Graph Clustering
        Raghavan et al. proposed a label propagation algorithm for graph clustering.
        This algorithm initializes every vertex of a graph with a unique label. Then, in iterative
        sweeps over the set of vertices the vertex labels are updated. A vertex gets the label
        that the maximum number of its neighbors have. The procedure is stopped when every vertex
        has the label that at least half of its neighbors have.
 *
 *
 */
class PLP: public NetworKit::CommunityDetectionAlgorithm {

protected:

	count updateThreshold = 0;
	count nIterations = 0; //!< number of iterations in last run


public:

	std::string VERSION;	// algorithm version number - increment in constructor for significant changes to the implementation

	PLP(count theta = none);

	virtual ~PLP();

	/**
	 * Run the label propagation clustering algorithm.
	 *
	 * @param[in]	G	input graph
	 * @return			clustering
	 */
	virtual Partition run(Graph& G);

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;


	/**
	 * The algorithm runs until a number of nodes less than
	 * the threshold is updated.
	 *
	 */
	virtual void setUpdateThreshold(count th);

	/**
	* Get number of iterations in last run.
	*/
	virtual count numberOfIterations();


};

} /* namespace NetworKit */
#endif /* PLP_H_ */
