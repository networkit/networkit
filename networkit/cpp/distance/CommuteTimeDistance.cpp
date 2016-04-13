/*
 * CommuteTimeDistance.cpp
 *
 *  Created on: 29.07.2015
 *      Author: henningm
 */


#include "CommuteTimeDistance.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Timer.h"

#include <fstream>
#include <sstream>
#include <math.h>

#include "omp.h"

namespace NetworKit {

CommuteTimeDistance::CommuteTimeDistance(const Graph& G, double tol): G(G), tol(tol), lamg(1e-5) {
	// prepare LAMG
	CSRMatrix matrix = CSRMatrix::graphLaplacian(G);
	Aux::Timer t;
	t.start();
	lamg.setupConnected(matrix);
	t.stop();

	setupTime = t.elapsedMilliseconds();

	DEBUG("done setting up Spanning");
}

void CommuteTimeDistance::run() {
	count n = G.numberOfNodes();
	distances.clear();
	distances.resize(n);
	G.forNodes([&](node v){
		distances[v].resize(n, 0.0);
	});

	// set up solution vector and status
	Vector solution(n);

	Vector rhs(n, 0.0);
	Vector zeroVector(n, 0.0);

	// solve for each pair of nodes
	G.forNodes([&](node u) {
		G.forNodes([&](node v) {
			// set up right-hand side
			rhs[u] = +1.0;
			rhs[v] = -1.0;
			TRACE("before solve for ", u, " and ", v);

			solution = zeroVector;

			lamg.solve(rhs, solution);
			double diff = solution[u] - solution[v];
			distances[u][v] = fabs(diff); // TODO: check unweighted, fix weighted case!
			rhs[u] = 0.0;
			rhs[v] = 0.0;
		});
	});
	hasRun = true;
}

uint64_t CommuteTimeDistance::getSetupTime() const {
	return setupTime;
}

void CommuteTimeDistance::runApproximation() {
	count n = G.numberOfNodes();
	distances.clear();
	distances.resize(n);
	G.forNodes([&](node v){
		distances[v].resize(n, 0.0);
	});
	double epsilon2 = tol * tol;
	const count k = ceil(log2(n)) / epsilon2;
	double randTab[3] = {1/sqrt(k), -1/sqrt(k)};
	Vector solution(n);

	for (index i = 0; i < k; ++i) {
		Vector rhs(n, 0.0);

		// matrix vector product of q
		// rhs(v) = \sum_e=1 ^m q(e) * B(e, v)
		//        = +/- q(e)
		G.forEdges([&](node u, node v) {
			double r = randTab[Aux::Random::integer(1)];

			if (u < v) {
				rhs[u] += r;
				rhs[v] -= r;
			}
			else {
				rhs[u] -= r;
				rhs[v] += r;
			}
		});

		lamg.solve(rhs, solution);

		G.forNodes([&](node u) {
			G.forNodes([&](node v) {
				double diff = solution[u] - solution[v];
				distances[u][v] += diff * diff; // TODO: fix weighted case!
			});
		});
	}

	hasRun = true;
}

void CommuteTimeDistance::runParallelApproximation() {
	count n = G.numberOfNodes();
	distances.clear();
	distances.resize(n);
	G.forNodes([&](node v){
		distances[v].resize(n, 0.0);
	});
	double epsilon2 = tol * tol;
	const count k = ceil(log2(n)) / epsilon2;
	double randTab[3] = {1/sqrt(k), -1/sqrt(k)};
	std::vector<Vector> solutions(k, Vector(n));
	std::vector<Vector> rhs(k, Vector(n));


#pragma omp parallel for
	for (index i = 0; i < k; ++i) {
		// rhs(v) = \sum_e=1 ^m q(e) * B(e, v)
		//        = +/- q(e)
		G.forEdges([&](node u, node v) {
			double r = randTab[Aux::Random::integer(1)];

			if (u < v) {
				rhs[i][u] += r;
				rhs[i][v] -= r;
			}
			else {
				rhs[i][u] -= r;
				rhs[i][v] += r;
			}
		});
	}

	lamg.parallelSolve(rhs, solutions);

	for (index i = 0; i < k; ++i) {
		G.parallelForNodes([&](node u) {
				G.forNodes([&](node v) {
				double diff = solutions[i][u] - solutions[i][v];
				distances[u][v] += diff * diff; // TODO: fix weighted case!
			});
		});
	}

	hasRun = true;
}

double CommuteTimeDistance::distance(node u, node v) {
	if (!hasRun) throw std::runtime_error("Call run method first");
	return sqrt(G.numberOfEdges()*distances[u][v]); // TODO fix weighted case: volume is the sum of the weights of the edges
}

}
