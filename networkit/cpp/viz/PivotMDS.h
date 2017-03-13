/*
 * PivotMDS.h
 *
 *  Created on: Jul 7, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_VIZ_PIVOTMDS_H_
#define NETWORKIT_CPP_VIZ_PIVOTMDS_H_

#include "GraphLayoutAlgorithm.h"

#include "../algebraic/CSRMatrix.h"
#include "../algebraic/Vector.h"

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup viz
 *
 * Implementation of PivotMDS proposed by Brandes and Pich.
 */
class PivotMDS : public GraphLayoutAlgorithm<double> {
public:
	/**
	 * Constructs a PivotMDS object for the given @a graph. The algorithm should embed the graph in @a dim dimensions
	 * using @a numPivots pivots.
	 * @param graph
	 * @param dim
	 * @param numPivots
	 */
	PivotMDS(const Graph& graph, count dim, count numPivots);

	/*
	 * Default destructor
	 */
	virtual ~PivotMDS() = default;

	/**
	 * Runs the PivotMDS algorithm.
	 */
	void run() override;

private:
	count dim;
	count numPivots;

	/**
	 *  Randomly picks the pivots for the algorithm
	 */
	std::vector<node> computePivots();

	/**
	 * Power method to compute the largest eigenvector and eigenvalue that are stored in @a eigenvector and
	 * @a eigenvalue.
	 */
	void powerMethod(const CSRMatrix& mat, const count n, Vector& eigenvector, double& eigenvalue);
};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_VIZ_PIVOTMDS_H_ */
