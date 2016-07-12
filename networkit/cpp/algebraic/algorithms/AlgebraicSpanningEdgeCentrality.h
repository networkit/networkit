/*
 * AlgebraicSpanningEdgeCentrality.h
 *
 *  Created on: Jul 12, 2016
 *      Author: Michael
 */

#ifndef NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICSPANNINGEDGECENTRALITY_H_
#define NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICSPANNINGEDGECENTRALITY_H_

#include "../../centrality/Centrality.h"
#include "../../numerics/LAMG/Lamg.h"

namespace NetworKit {

/**
 * @ingroup algebraic
 * Implementation of Spanning edge centrality with algebraic notation.
 */
template<class Matrix>
class AlgebraicSpanningEdgeCentrality : public Centrality {
public:
	/**
	 * Constructs an instance of the AlgebraicSpanningEdgeCentrality algorithm for the given Graph @a graph.
	 * The tolerance @a tol is used to control the approximation error when approximating the spanning edge
	 * centrality for @a graph.
	 * @param graph
	 * @param tol
	 */
	AlgebraicSpanningEdgeCentrality(const Graph& graph, double tol = 0.1) : Centrality(graph), tol(tol) {}
	virtual ~AlgebraicSpanningEdgeCentrality() = default;


	/**
	 * Compute spanning edge centrality exactly.
	 */
	void run() override;

	/**
	 * Approximate spanning edge centrality scores with the Johnson-Lindenstrauss transform.
	 */
	void runApproximation();

private:
	double tol;
};

template<class Matrix>
void AlgebraicSpanningEdgeCentrality<Matrix>::run() {
	const count n = G.numberOfNodes();
	const count m = G.numberOfEdges();
	scoreData.clear();
	scoreData.resize(m, 0.0);


	std::vector<Vector> rhs(m, Vector(n));
	this->G.parallelForEdges([&](node u, node v, edgeid e) {
		rhs[e][u] = +1;
		rhs[e][v] = -1;
	});

	std::vector<Vector> solutions(m, Vector(n));

	Lamg<Matrix> lamg(1e-5);
	lamg.setupConnected(Matrix::laplacianMatrix(this->G));
	lamg.parallelSolve(rhs, solutions);

	this->G.parallelForEdges([&](node u, node v, edgeid e) {
		double diff = solutions[e][u] - solutions[e][v];
		scoreData[e] = fabs(diff);
	});

	hasRun = true;
}

template<class Matrix>
void AlgebraicSpanningEdgeCentrality<Matrix>::runApproximation() {
	const count n = G.numberOfNodes();
	const count m = G.numberOfEdges();
	scoreData.clear();
	scoreData.resize(m, 0.0);
	double epsilon2 = tol * tol;
	const count k = ceil(log2(n)) / epsilon2;
	double randTab[3] = {1.0/sqrt(k), -1.0/sqrt(k)};

	// compute random projection matrix Q
	std::vector<Triplet> triplets(k*m);
	index tripletIdx = 0;
	for (index i = 0; i < k; ++i) {
		for (index j = 0; j < m; ++j) {
			double val = randTab[Aux::Random::integer(1)];
			triplets[tripletIdx++] = {i,j,val};
		}
	}

	Matrix Q(k, m, triplets);
	Matrix B = Matrix::incidenceMatrix(this->G).transpose(); // transposed incidence matrix with dimensions m x n

	Matrix Y = Q * B; // k x n

	std::vector<Vector> yRows(k, Vector(n));
	std::vector<Vector> zRows(k, Vector(n));
#pragma omp parallel for
	for (index i = 0; i < k; ++i) {
		yRows[i] = Y.row(i);
	}

	Lamg<Matrix> lamg(1e-5);
	Matrix laplacian = Matrix::laplacianMatrix(this->G);
	lamg.setupConnected(laplacian);
	lamg.parallelSolve(yRows, zRows);

	std::vector<Vector> zColumns(n, Vector(k));
#pragma omp parallel for
	for (index j = 0; j < n; ++j) {
		for (index i = 0; i < k; ++i) {
			zColumns[j][i] = zRows[i][j];
		}
	}

	this->G.parallelForEdges([&](node u, node v, edgeid e) {
		double norm = (zColumns[u] - zColumns[v]).length();
		scoreData[e] = norm * norm;
	});

	hasRun = true;
}


} /* namespace NetworKit */



#endif /* NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICSPANNINGEDGECENTRALITY_H_ */
