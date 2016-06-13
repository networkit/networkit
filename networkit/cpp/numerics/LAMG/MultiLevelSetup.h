/*
 * MultiLevelSetup.h
 *
 *  Created on: 10.01.2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef MULTILEVELSETUP_H_
#define MULTILEVELSETUP_H_

#include "LevelHierarchy.h"
#include "../Smoother.h"
#include "../../algebraic/CSRMatrix.h"

#include <limits>

namespace NetworKit {

#define UNDECIDED std::numeric_limits<index>::max()

/**
 * @ingroup
 * Implements the setup phase of LAMG (Lean Algebraic Multigrid by Livne et al.).
 */
class MultiLevelSetup {

private:
	const Smoother &smoother;
#ifndef NDEBUG
	static count eliminationTime;
	static count schurComplementTime;
	static count aggregationTime;
#endif

	/**
	 * Elimination phase of LAMG for the specified Laplacian matrix @a matrix. If an elemination stage is created,
	 * it is stored in the LevelHierarchy @a hierarchy and the method returns @code{true}. Otherwise, @code{false}
	 * is returned.
	 * @param matrix Laplacian matrix for which an elimination stage should be created.
	 * @param hierarchy LevelHierarchy in which the created elimination stage is stored (if created)
	 * @return @code{True} if an elimination stage has been created, otherwise @code{false}.
	 */
	bool coarseningElimination(CSRMatrix &matrix, LevelHierarchy &hierarchy) const;

	/**
	 * Scans the Laplacian matrix for nodes with a low degree (i.e. nodes with less than 5 neighbors). For each node
	 * having low degree and independent from already found low degree nodes, @code{true} is stored in @a fNode. The
	 * @a stage parameter specifies if we are in the first or subsequent stages during elimination.
	 * @param matrix Laplacian matrix.
	 * @param fNode[out] For each node, @code{true} if the node is of low degree and @code{false} otherwise.
	 * @param stage The stage of the elimination phase.
	 * @return Number of low degree nodes found.
	 */
	count lowDegreeSweep(const CSRMatrix &matrix, std::vector<bool> &fNode, index stage) const;

	/**
	 * Computes the projection matrix @a P and the @a q vector used to restrict and interpolate the matrix for an
	 * elimination stage.
	 * @param matrix Laplacian matrix.
	 * @param fSet Set of nodes having low degree.
	 * @param coarseIndex Set of nodes equal to V \setminus fSet
	 * @param P[out] The projection matrix.
	 * @param q[out] The q vector.
	 */
	void eliminationOperators(const CSRMatrix &matrix, const std::vector<index> &fSet, const std::vector<index> &coarseIndex, CSRMatrix &P, Vector &q) const;

	/**
	 * Computes the submatrix specified by @a rows and @a columns and stores the result in @a result. The @a coarseIndex
	 * is used to efficiently check whether the value is part of the submatrix or not.
	 * @param matrix Laplacian matrix.
	 * @param rows Rows indices to include in submatrix.
	 * @param columns Column indices to include in submatrix.
	 * @param coarseIndex
	 * @param result[out] The resulting submatrix.
	 */
	void subMatrix(const CSRMatrix &matrix, const std::vector<index> &rows, const std::vector<index> &columns, const std::vector<index> &coarseIndex, CSRMatrix &result) const;


	/**
	 * Aggregation phase of LAMG for the specified Laplacian matrix @a matrix. The coarsened matrix is stored in the
	 * LevelHierarchy @a hierarchy. The test vector @a tv is used for determining the affinity of nodes. The parameter
	 * @a numTVVectors specifies how many test vectors will be created to determine the affinities between nodes.
	 * @param matrix Laplacian matrix.
	 * @param hierarchy LevelHierarchy to store the stage.
	 * @param tv Test vector.
	 * @param numTVVectors Number of test vectors to use.
	 */
	void coarseningAggregation(CSRMatrix &matrix, LevelHierarchy &hierarchy, Vector &tv, count numTVVectors) const;

	/**
	 * Create @a numVectors test vectors for the given Laplacian matrix @a matrix. The test vector @a tv will be
	 * reused.
	 * @param matrix Laplacian matrix.
	 * @param tv Test vector.
	 * @param numVectors Number of test vectors to create.
	 * @return The created test vectors.
	 */
	std::vector<Vector> generateTVs(const CSRMatrix &matrix, Vector &tv, const count numVectors) const;

	/**
	 * Adds high degree nodes as seeds to @a status.
	 * @param matrix Laplacian matrix.
	 * @param status[out] High degree nodes are flagged as seed.
	 */
	void addHighDegreeSeedNodes(const CSRMatrix &matrix, std::vector<index> &status) const;

	/**
	 * Aggregates loose nodes (nodes with low adjacency) together.
	 * @param strongAdjMatrix Strong adjacency matrix.
	 * @param status[out] Aggregates loose nodes together and labels them in @a status accordingly.
	 * @param nc[out] The altered number of coarse nodes.
	 */
	void aggregateLooseNodes(const CSRMatrix &strongAdjMatrix, std::vector<index> &status, count &nc) const;

	/**
	 * Computes the strong adjacency matrix for the given Laplacian matrix @a matrix.
	 * @param matrix Laplacian matrix.
	 * @param strongAdjMatrix[out] The resulting strong adjacency matrix.
	 */
	void computeStrongAdjacencyMatrix(const CSRMatrix &matrix, CSRMatrix &strongAdjMatrix) const;

	/**
	 * Computes the affinity matrix for the given Laplacian matrix @a matrix and the test vectors @a tVs.
	 * @param matrix Laplacian matrix.
	 * @param tVs Test vectors.
	 * @param affinityMatrix[out] The resulting affinity matrix.
	 */
	void computeAffinityMatrix(const CSRMatrix &matrix, const std::vector<Vector> &tVs, CSRMatrix &affinityMatrix) const;

	/**
	 * Models one stage in the aggregation phase. New aggregates are labeled accordingly in @a status.
	 * @param matrix Laplacian matrix.
	 * @param nc Number of coarse nodes.
	 * @param strongAdjMatrix Strong adjacency matrix.
	 * @param affinityMatrix Affinity matrix.
	 * @param tVs[out] Test vectors.
	 * @param status[out] Aggregation labels.
	 */
	void aggregationStage(const CSRMatrix &matrix, count &nc, const CSRMatrix &strongAdjMatrix, const CSRMatrix &affinityMatrix, std::vector<Vector> &tVs, std::vector<index> &status) const;

	/**
	 * Computes strong (cf. LAMG paper) neighbors and stores them sorted into @a bins.
	 * @param affinityMatrix Affinity matrix.
	 * @param status Aggregation labels.
	 * @param bins[out] Strong neighbors sorted into bins.
	 */
	void computeStrongNeighbors(const CSRMatrix &affinityMatrix, const std::vector<index> &status, std::vector<std::vector<index>> &bins) const;

	/**
	 * Finds the best seed for node @a u and stores it in @a s.
	 * @param strongAdjMatrix
	 * @param affinityMatrix Affinity matrix.
	 * @param diag Vector of diagonal entries of the Laplacian matrix.
	 * @param tVs Test vectors.
	 * @param status Aggregation labels.
	 * @param u The node to find the best seed for.
	 * @param s[out] The best seed for node @a u.
	 * @return @code{True} if a seed has been found for @a u, @code{false} otherwise.
	 */
	bool findBestSeedEnergyCorrected(const CSRMatrix &strongAdjMatrix, const CSRMatrix &affinityMatrix, const std::vector<double> &diag, const std::vector<Vector> &tVs, const std::vector<index> &status, const index u, index &s) const;


	/**
	 * Determines if the Laplacian matrix @a A can be coarsened further.
	 * @param A Laplacian matrix.
	 * @return @code{True} if @a A can be coarsened further, @code{false} otherwise.
	 */
	bool canCoarsen(const CSRMatrix &A) const;

	/**
	 * Determines if the relaxation is fast enough to stop coarsening.
	 * @param A Laplacian matrix.
	 * @param lvlIndex The number of levels already created in the hierarchy.
	 * @param tv Test vector.
	 * @return @code{True} if convergence of relaxation is fast, @code{false} otherwise.
	 */
	bool isRelaxationFast(const CSRMatrix &A, index lvlIndex, Vector &tv) const;

	/**
	 * Computes the coarsened matrix of @a matrix by means of the projection matrix @a P and stores the result in @a B.
	 * @param P Projection matrix.
	 * @param A Laplacian matrix.
	 * @param PColIndex Stores the column index of the 1 entry at each row.
	 * @param PRowIndex Stores the row index of the 1 entry at each column.
	 * @param B[out] Resulting coarsened Laplacian matrix.
	 */
	void galerkinOperator(const CSRMatrix &P, const CSRMatrix &A, const std::vector<index> &PColIndex, const std::vector<std::vector<index>> &PRowIndex, CSRMatrix &B) const;

public:
	/**
	 * Creates an instance of MultiLevelSetup with the specified @a smoother used for relaxing during the setup phase.
	 * @param smoother Reference to smoother.
	 */
	MultiLevelSetup(const Smoother &smoother);

	/**
	 * Creates a @å hierarchy for the given Laplacian matrix of the graph @a G.
	 * @param G The graph.
	 * @param hierarchy[out] The constructed hierarchy.
	 */
	void setup(const Graph &G, LevelHierarchy &hierarchy) const;

	/**
	 * Creates a @a hierarchy for the given Laplacian matrix @a matrix.
	 * @param matrix Laplcian matrix.
	 * @param hierarchy[out] The constructed hierarchy.
	 */
	void setup(const CSRMatrix &matrix, LevelHierarchy &hierarchy) const;
};

} /* namespace NetworKit */

#endif /* MULTILEVELSETUP_H_ */
