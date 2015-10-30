/*
 * Spanning.cpp
 *
 *  Created on: 29.07.2015
 *      Author: henningm
 */


#include "Spanning.h"
#include "../auxiliary/Log.h"
#include "../numerics/LAMG/SolverLamg.h"
#include "../numerics/LAMG/MultiLevelSetup.h"
#include "../numerics/GaussSeidelRelaxation.h"
#include "../algebraic/LaplacianMatrix.h"
#include "../spanning/RandomSpanningTree.h"
#include "../spanning/PseudoRandomSpanningTree.h"

namespace NetworKit {

Spanning::Spanning(const Graph& G, double tol): Centrality(G), tol(tol) {
	// prepare LAMG
	smoother = new GaussSeidelRelaxation(tol);
	CSRMatrix matrix = CSRMatrix::graphLaplacian(G);

	MultiLevelSetup setup(*smoother);
	hierarchy = new LevelHierarchy();
	setup.setup(matrix, *hierarchy);
	solver = new SolverLamg(*hierarchy, *smoother);
	DEBUG("done setting up Spanning");
}

Spanning::~Spanning() {
	delete smoother;
	delete solver;
	delete hierarchy;
}


void Spanning::run() {
	count n = G.numberOfNodes();
	scoreData.clear();
	scoreData.resize(G.numberOfEdges(), 0.0);

	// set up solution vector and status
	Vector solution(n);
	LAMGSolverStatus status;
	status.desiredResidualReduction = tol;

	Vector rhs(n, 0.0);
	Vector zeroVector(n, 0.0);

	// solve for each edge
	G.forEdges([&](node u, node v, edgeid e) {
		// set up right-hand side
		rhs[u] = +1.0;
		rhs[v] = -1.0;
		TRACE("before solve for ", u, " and ", v);

		solution = zeroVector;

		solver->solve(solution, rhs, status);
		double diff = solution[u] - solution[v];
		scoreData[e] = fabs(diff); // TODO: check unweighted, fix weighted case, fix edge IDs!
		rhs[u] = 0.0;
		rhs[v] = 0.0;
	});
}


void Spanning::runApproximation() {
	const count n = G.numberOfNodes();
	const count m = G.numberOfEdges();
	const count k = 8 * ceil(log(n));
	double randTab[3] = {0, 1/sqrt(k), -1/sqrt(k)};
	LAMGSolverStatus status;
	status.desiredResidualReduction = tol;
	Vector solution(n);
	scoreData.clear();
	scoreData.resize(m, 0.0);

	for (index i = 0; i < k; ++i) {
		Vector rhs(n, 0.0);

		// rhs(v) = \sum_e=1 ^m q(e) * B(e, v)
		//        = +/- q(e)
		G.forEdges([&](node u, node v) {
			double r = randTab[Aux::Random::integer(2)];
			if (u < v) {
				rhs[u] += r;
				rhs[v] -= r;
			}
			else {
				rhs[u] -= r;
				rhs[v] += r;
			}
		});

		solver->solve(solution, rhs, status);

		G.forEdges([&](node u, node v, edgeid e) {
			double diff = solution[u] - solution[v];
			scoreData[e] += diff * diff; // TODO: fix weighted case!
		});
	}
}

void Spanning::runTreeApproximation() {
	const count reps = 500;
	scoreData.clear();
	scoreData.resize(G.numberOfEdges(), 0.0);

	RandomSpanningTree rst(G);

	for (index i = 0; i < reps; ++i) {
		rst.run();
		Graph tree = rst.getTree();
		G.forEdges([&](node u, node v, edgeid e) {
			scoreData[e] += tree.hasEdge(u, v);
		});
	}

	G.forEdges([&](node u, node v, edgeid e) {
		scoreData[e] /= reps;
	});
}

void Spanning::runPseudoTreeApproximation() {
	const count reps = 500;
	scoreData.clear();
	scoreData.resize(G.numberOfEdges(), 0.0);

	PseudoRandomSpanningTree rst(G);

	for (index i = 0; i < reps; ++i) {
		rst.run();
		Graph tree = rst.getTree();
		G.forEdges([&](node u, node v, edgeid e) {
			scoreData[e] += tree.hasEdge(u, v);
		});
	}

	G.forEdges([&](node u, node v, edgeid e) {
		scoreData[e] /= reps;
	});
}

double Spanning::runForEdge(node u, node v) {
	count n = G.numberOfNodes();

	// set up solution vector and status
	Vector solution(n, 0.0);
	Vector rhs(n, 0.0);

	LAMGSolverStatus status;
	status.desiredResidualReduction = tol;

	// set up right-hand side
	rhs[u] = +1.0;
	rhs[v] = -1.0;
	TRACE("before solve for ", u, " and ", v);

	solver->solve(solution, rhs, status);
	return fabs(solution[u] - solution[v]); // TODO: fix weighted case, fix edge IDs!
}


}
