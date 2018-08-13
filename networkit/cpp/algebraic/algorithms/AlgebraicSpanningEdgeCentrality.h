/*
 * AlgebraicSpanningEdgeCentrality.h
 *
 *  Created on: Jul 12, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
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

	std::vector<Vector> yRows(k, Vector(n));
	std::vector<Vector> zRows(k, Vector(n));

	G.forEdges([&](node u, node v, edgeweight w, index) {
#pragma omp parallel for
		for (omp_index i = 0; i < static_cast<omp_index>(k); ++i) {
			double rand = randTab[Aux::Random::integer(1)];
			yRows[i][u] += w * rand;
			yRows[i][v] -= w * rand;
		}
	});


	Lamg<Matrix> lamg(1e-5);
	lamg.setupConnected(Matrix::laplacianMatrix(this->G));
	lamg.parallelSolve(yRows, zRows);

	this->G.parallelForEdges([&](node u, node v, edgeid e) {
		double sqSum = 0.0;
		for (index i = 0; i < k; ++i) {
			double diff = (zRows[i][u] - zRows[i][v]);
			sqSum += diff * diff;
		}

		scoreData[e] = sqSum;
	});

	hasRun = true;
}


} /* namespace NetworKit */



#endif /* NETWORKIT_CPP_ALGEBRAIC_ALGORITHMS_ALGEBRAICSPANNINGEDGECENTRALITY_H_ */
