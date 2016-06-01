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

CommuteTimeDistance::CommuteTimeDistance(const Graph& G, double tol): Algorithm(), G(G), tol(tol), lamg(1e-5) {
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
	G.forNodePairs([&](node u, node v){
		// set up right-hand side
		rhs[u] = +1.0;
		rhs[v] = -1.0;
		TRACE("before solve for ", u, " and ", v);

		solution = zeroVector;

		lamg.solve(rhs, solution);
		double diff = solution[u] - solution[v];
		distances[u][v] = fabs(diff); // TODO: check unweighted, fix weighted case!
		distances[v][u] = fabs(diff); // TODO: check unweighted, fix weighted case!
		rhs[u] = 0.0;
		rhs[v] = 0.0;
	});
	exactly = true;
	hasRun = true;
}

uint64_t CommuteTimeDistance::getSetupTime() const {
	return setupTime;
}

void CommuteTimeDistance::runApproximation() {
	count n = G.numberOfNodes();
	// distances.clear();
	// distances.resize(n);
	// G.forNodes([&](node v){
	// 	distances[v].resize(n, 0.0);
	// });
	double epsilon2 = tol * tol;
	k = ceil(log2(n)) / epsilon2;
	double randTab[3] = {1/sqrt(k), -1/sqrt(k)};
	solutions.clear();
	solutions.resize(k, Vector(n));

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

		lamg.solve(rhs, solutions[i]);

		// G.forNodePairs([&](node u, node v){
		// 		double diff = solutions[i][u] - solutions[i][v];
		// 		distances[u][v] += diff * diff; // TODO: fix weighted case!
		// 		distances[v][u] += diff * diff; // TODO: fix weighted case!
		// });
	}
	exactly = false;
	hasRun = true;
}

void CommuteTimeDistance::runParallelApproximation() {
	count n = G.numberOfNodes();
	// distances.clear();
	// distances.resize(n);
	// G.forNodes([&](node v){
	// 	distances[v].resize(n, 0.0);
	// });
	double epsilon2 = tol * tol;
	k = ceil(log2(n)) / epsilon2;
	double randTab[3] = {1/sqrt(k), -1/sqrt(k)};
	solutions.clear();
	solutions.resize(k, Vector(n));
	std::vector<Vector> rhs(k, Vector(n));

	INFO("Number k of iterations: ", k);
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
	INFO("Starting the solve phase");
	lamg.parallelSolve(rhs, solutions);
	INFO("Done with the solve phase");

	// for (index i = 0; i < k; ++i) {
	// 	G.parallelForNodePairs([&](node u, node v){
	// 		double diff = solutions[i][u] - solutions[i][v];
	// 		distances[u][v] += diff * diff; // TODO: fix weighted case!
	// 		distances[v][u] += diff * diff; // TODO: fix weighted case!
	// 	});
	// }
	exactly = false;
	hasRun = true;
}

double CommuteTimeDistance::distance(node u, node v) {
	if (!hasRun) throw std::runtime_error("Call run method first");
	if (exactly) {
		return sqrt(distances[u][v]* G.numberOfEdges()); // TODO fix weighted case: volume is the sum of the weights of the edges
	} else {
		double dist = 0;
		for (index i = 0; i < k; ++i) {
			double diff = solutions[i][u] - solutions[i][v];
			dist += diff * diff;
		}
		return sqrt(dist* G.numberOfEdges());
	}
}

double CommuteTimeDistance::runSinglePair(node u, node v) {
	count n = G.numberOfNodes();
	double dist = 0.0;

	// set up solution vector and status
	Vector solution(n);

	Vector rhs(n, 0.0);
	Vector zeroVector(n, 0.0);
	rhs[u] = +1.0;
	rhs[v] = -1.0;
	// set up right-hand side
	solution = zeroVector;
	lamg.solve(rhs, solution);
	double diff = solution[u] - solution[v];
	dist = fabs(diff); // TODO: check unweighted, fix weighted case!
	return sqrt(dist* G.numberOfEdges());
}

double CommuteTimeDistance::runSingleSource(node u) {
	count n = G.numberOfNodes();
	double dist = 0.0;
	double sum = 0.0;
	Vector zeroVector(n, 0.0);
	// set up solution vector and status
	std::vector<Vector> rhs(n, Vector(n));
	std::vector<Vector> solution(n, Vector(n));
	G.forNodes([&](node i){
		rhs[i] = zeroVector;
		solution[i] = zeroVector;
		rhs[i][u] = +1.0;
		if (i != u) {
			rhs[i][i] = -1.0;
		} else {
			rhs[i][0] = -1.0;
		}
	});
	INFO("rhs.size() = ", rhs.size());
	INFO("solutions.size() = ", solution.size());
	lamg.parallelSolve(rhs, solution);
	G.forNodes([&](node i){
		if (i != u) {
			double diff = solution[i][u] - solution[i][i];
			dist = fabs(diff); // TODO: check unweighted, fix weighted case!
			sum += sqrt(dist);
		}
	});
	return sum * sqrt(G.numberOfEdges());
}

}
