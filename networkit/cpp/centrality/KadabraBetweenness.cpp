/*
 * KadabraBetweenness.cpp
 *
 * Created on: 18.07.2018
 * 		 Author: Eugenio Angriman, Alexander van der Grinten
 */

#include <cmath>
#include <ctime>
#include <omp.h>

#include "../auxiliary/Random.h"
#include "../auxiliary/Timer.h"
#include "../distance/Diameter.h"
#include "KadabraBetweenness.h"

namespace NetworKit {

Status::Status(const count k)
    : k(k), top(k), approxTop(k), finished(k), bet(k), errL(k), errU(k) {}

KadabraBetweenness::KadabraBetweenness(const Graph &G, const double err,
                                       const double delta, const count k,
                                       count unionSample,
                                       const count startFactor)
    : G(G), err(err), delta(delta), k(k), n(G.upperNodeIdBound()),
      startFactor(startFactor), unionSample(unionSample), absolute(k == 0) {
	if (k > n) {
		throw std::runtime_error(
		    "k is higher than the number of nodes of the input graph! Choose a "
		    "value between 0 (absolute) and n-1.");
	}

	if (delta >= 1 || delta <= 0) {
		throw std::runtime_error(
		    "Delta should be greater than 0 and smaller than 1.");
	}

	if (err >= 1 || err <= 0) {
		throw std::runtime_error(
		    "The error should be greater than 0 and smaller than 1.");
	}

	if (!Aux::Random::getUseThreadId()) {
		throw std::runtime_error(
		    "Error: the Kadabra algorithm needs 'useThreadId' set to true.");
	}
}

bool KadabraBetweenness::computeFinished(Status *status) const {
	std::vector<double> &bet = status->bet;
	std::vector<double> &errL = status->errL;
	std::vector<double> &errU = status->errU;
	bool allFinished = true;

	count i;
	for (i = 0; i < status->k - 1; ++i) {
		bet[i] = status->approxTop[i] / (double)status->nPairs;
		errL[i] = computeF(bet[i], status->nPairs, deltaLGuess[status->top[i]]);
		errU[i] = computeG(bet[i], status->nPairs, deltaUGuess[status->top[i]]);
	}

	bet[i] = status->approxTop[i] / (double)status->nPairs;
	errL[i] = computeF(bet[i], status->nPairs, this->deltaLMinGuess);
	errU[i] = computeG(bet[i], status->nPairs, this->deltaUMinGuess);

	if (absolute) {
		for (count i = 0; i < status->k; ++i) {
			status->finished[i] = (errL[i] < err && errU[i] < err);
			allFinished = allFinished && status->finished[i];
		}
	} else {
		for (count i = 0; i < status->k; ++i) {
			if (i == 0) {
				status->finished[i] = (bet[i] - errL[i] > bet[i + 1] + errU[i + 1]);
			} else if (i < k) {
				status->finished[i] = (bet[i - 1] - errL[i - 1] > bet[i] + errU[i]) &&
				                      (bet[i] - errL[i] > bet[i + 1] + errU[i + 1]);
			} else {
				status->finished[i] = bet[k - 1] - errU[k - 1] > bet[i] + errU[i];
			}
			status->finished[i] =
			    status->finished[i] || (errL[i] < err && errU[i] < err);
			allFinished = allFinished && status->finished[i];
		}
	}

	return allFinished;
}

// Computes the function f that bounds the betweenness of a vertex from below.
// For more information, see Borassi, Natale (2016).
double KadabraBetweenness::computeF(const double btilde, const count iterNum,
                                    const double deltaL) const {
	double tmp = (((double)omega) / iterNum - 1. / 3);
	double errChern = (std::log(1. / deltaL)) * 1. / iterNum *
	                  (-tmp + std::sqrt(tmp * tmp + 2 * btilde * omega /
	                                                    (std::log(1. / deltaL))));
	return std::min(errChern, btilde);
}

// Computes the function g that bounds the betweenness of a vertex from above.
// For more information, see Borassi, Natale (2016).
double KadabraBetweenness::computeG(const double btilde, const count iterNum,
                                    const double deltaU) const {
	double tmp = (((double)omega) / iterNum + 1. / 3);
	double errChern = (std::log(1. / deltaU)) * 1. / iterNum *
	                  (tmp + std::sqrt(tmp * tmp + 2 * btilde * omega /
	                                                   (std::log(1. / deltaU))));
	return std::min(errChern, 1 - btilde);
}

void KadabraBetweenness::oneRound(SpSampler &sampler) {
	auto path = sampler.randomPath();
	for (node u : path) {
		approx[omp_get_thread_num()][u] += 1.;
	}
}

void KadabraBetweenness::getStatus(Status *status, const bool parallel) const {
	if (status != NULL) {
		auto loop = [&](count i) {
			if (absolute) {
				status->top[i] = i;
				status->approxTop[i] = approxSum[i];
			} else {
				status->top[i] = top->getElement(i);
				status->approxTop[i] = top->getValue(i);
			}
		};
		if (parallel) {
#pragma omp parallel for
			for (omp_index i = 0; i < static_cast<omp_index>(unionSample); ++i) {
				loop(static_cast<count>(i));
			}
		} else {
			for (count i = 0; i < unionSample; ++i) {
				loop(i);
			}
		}
		status->nPairs = nPairs;
	}
}

void KadabraBetweenness::computeBetErr(Status *status, std::vector<double> &bet,
                                       std::vector<double> &errL,
                                       std::vector<double> &errU) const {
	count i;
	double maxErr = std::sqrt(startFactor) * err / 4.;

	for (i = 0; i < status->k; ++i) {
		bet[i] = status->approxTop[i] / (double)status->nPairs;
	}

	if (absolute) {
		for (i = 0; i < status->k; ++i) {
			errL[i] = err;
			errU[i] = err;
		}
	} else {
		errU[0] = std::max(err, (bet[0] - bet[1]) / 2.);
		errL[0] = 10;
		for (i = 1; i < k; ++i) {
			errL[i] = std::max(err, (bet[i - 1] - bet[i]) / 2.);
			errU[i] = std::max(err, (bet[i] - bet[i + 1]) / 2.);
		}
		for (i = k; i < status->k; ++i) {
			errL[i] = 10;
			errU[i] = std::max(err, bet[k - 1] + (bet[k - 1] - bet[k]) / 2. - bet[i]);
		}
		for (i = 0; i < k - 1; ++i) {
			if (bet[i] - bet[i + 1] < maxErr) {
				errL[i] = err;
				errU[i] = err;
				errL[i + 1] = err;
				errU[i + 1] = err;
			}
		}
		for (i = k + 1; i < status->k; ++i) {
			if (bet[k] - bet[i] < maxErr) {
				errL[k] = err;
				errU[k] = err;
				errL[i] = err;
				errU[i] = err;
			}
		}
	}
}

void KadabraBetweenness::computeDeltaGuess() {
	double a = 0,
	       b = 1. / err / err * std::log(n * 4 * (1 - balancingFactor) / delta),
	       c = (a + b) / 2;
	double sum;

	Status status(unionSample);
	getStatus(&status, true);

	std::vector<double> bet(status.k);
	std::vector<double> errL(status.k);
	std::vector<double> errU(status.k);

	computeBetErr(&status, bet, errL, errU);

	for (count i = 0; i < unionSample; ++i) {
		count v = status.top[i];
		approxSum[v] = approxSum[v] / (double)nPairs;
	}

	while (b - a > err / 10.) {
		c = (b + a) / 2.;
		sum = 0;
		for (count i = 0; i < unionSample; ++i) {
			sum += std::exp(-c * errL[i] * errL[i] / bet[i]);
			sum += std::exp(-c * errU[i] * errU[i] / bet[i]);
		}

		sum += std::exp(-c * errL[unionSample - 1] * errL[unionSample - 1] /
		                bet[unionSample - 1]) *
		       (n - unionSample);
		sum += std::exp(-c * errU[unionSample - 1] * errU[unionSample - 1] /
		                bet[unionSample - 1]) *
		       (n - unionSample);

		if (sum >= delta / 2. * (1 - balancingFactor)) {
			a = c;
		} else {
			b = c;
		}
	}

	deltaLMinGuess = std::exp(-b * errL[unionSample - 1] * errL[unionSample - 1] /
	                          bet[unionSample - 1]) +
	                 delta * balancingFactor / 4. / (double)n;
	deltaUMinGuess = std::exp(-b * errU[unionSample - 1] * errU[unionSample - 1] /
	                          bet[unionSample - 1]) +
	                 delta * balancingFactor / 4. / (double)n;

	std::fill(deltaLGuess.begin(), deltaLGuess.end(), deltaLMinGuess);
	std::fill(deltaUGuess.begin(), deltaUGuess.end(), deltaUMinGuess);

	for (count i = 0; i < unionSample; ++i) {
		node v = status.top[i];
		deltaLGuess[v] = std::exp(-b * errL[i] * errL[i] / bet[i]) +
		                 delta * balancingFactor / 4. / (double)n;
		deltaUGuess[v] = std::exp(-b * errU[i] * errU[i] / bet[i]) +
		                 delta * balancingFactor / 4. / (double)n;
	}
}

void KadabraBetweenness::computeApproxParallel(const bool normalize) {
#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
		approxSum[i] = 0.;
		for (count j = 0; j < omp_max_threads; ++j) {
			approxSum[i] += approx[j][i];
		}
		if (normalize) {
			approxSum[i] /= (double)nPairs;
			if (!G.isDirected()) {
				approxSum[i] *= 2.;
			}
		}
	}
}

void KadabraBetweenness::init() {
	omp_max_threads = omp_get_max_threads();
	approx.assign(omp_max_threads, std::vector<double>(n, 0.));
	approxSum.resize(n);
	deltaLGuess.resize(n);
	deltaUGuess.resize(n);
	nPairs = 0;
	if (!G.isDirected()) {
		cc = new ConnectedComponents(G);
		cc->run();
	}
}

void KadabraBetweenness::fillResult() {
	if (absolute) {
		topkScores.resize(n);
		topkNodes.resize(n);
		rankingVector.resize(n);
#pragma omp parallel for
		for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
			rankingVector[i] = std::make_pair(i, approxSum[i]);
		}
		std::sort(rankingVector.begin(), rankingVector.end(),
		          [&](std::pair<node, double> p1, std::pair<node, double> p2) {
			          return p1.second > p2.second;
		          });
#pragma omp parallel for
		for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
			topkNodes[i] = rankingVector[i].first;
			topkScores[i] = rankingVector[i].second;
		}
	} else {
		topkScores.resize(k);
		topkNodes.resize(k);
		rankingVector.resize(k);
		for (count i = 0; i < k; ++i) {
			topkNodes[i] = top->getElement(i);
			topkScores[i] = approxSum[topkNodes[i]];
			assert(top->getValue(i) == topkScores[i]);
			rankingVector[i] = std::make_pair(topkNodes[i], topkScores[i]);
		}
	}
}

void KadabraBetweenness::run() {
	init();

	// TODO: setting the maximum relateve error to 0 gives the exact diameter but
	// may be inefficient for large graphs. What is the maximum relative error
	// that we can tolerate?
	Diameter diam(G, estimatedRange, 0.f);
	diam.run();
	// Getting diameter upper bound
	int32_t diameter = diam.getDiameter().second;
	omega =
	    0.5 / err / err * (std::log2(diameter - 1) + 1 + std::log(0.5 / delta));

	const count tau = omega / startFactor;

	if (unionSample == 0) {
		// In the absolute case we need to check that all the estimated betweenness
		// scores are within the error bounds. Thus, we set unionSample to the
		// number of nodes.
		if (absolute) {
			unionSample = n;
		} else {
			unionSample = std::min(
			    n,
			    (count)std::max((2 * std::sqrt(G.numberOfEdges()) / omp_max_threads),
			                    k + 20.));
		}
	}

	if (!absolute) {
		this->top = new Aux::SortedList(unionSample, n);
	}

#pragma omp parallel
	{
		SpSampler sampler(G, *cc);
		while (nPairs <= tau) {
			oneRound(sampler);
			++nPairs;
		}
	}

	computeApproxParallel();
	if (!absolute) {
		fillPQ();
	}
	computeDeltaGuess();
	nPairs = 0;
	std::atomic<bool> stop(false);
	if (!absolute) {
		top->clear();
	}
#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(omp_max_threads); ++i) {
		std::fill(approx[i].begin(), approx[i].end(), 0.);
	}

#pragma omp parallel
	{
		SpSampler sampler(G, *cc);
		Status status(unionSample);
		status.nPairs = 0;

		while (!stop && nPairs < omega) {
			for (unsigned short i = 0; i < itersPerStep; ++i) {
				oneRound(sampler);
			}
			nPairs += itersPerStep;
			if (omp_get_thread_num() == 0) {
				for (count i = 0; i < n; ++i) {
					approxSum[i] = 0.;
					for (count j = 0; j < omp_max_threads; ++j) {
						approxSum[i] += approx[j][i];
					}
					if (!absolute) {
						top->insert(i, approxSum[i]);
					}
				}

				getStatus(&status);
				stop = computeFinished(&status);
			}
		}
	}

	computeApproxParallel(true);
	if (!absolute) {
		// It should not be necessary to clear it again, but otherwise the
		// ranking is wrong.
		top->clear();
		fillPQ();
	}
	fillResult();
	nPairs += tau;
	hasRun = true;
	if (!absolute) {
		delete (top);
	}
}

SpSampler::SpSampler(const Graph &G, const ConnectedComponents &cc)
    : G(G), n(G.upperNodeIdBound()), pred(n, false, true), cc(cc) {
	q.resize(n);
	ballInd.assign(n, 0);
	dist.resize(n);
	nPaths.resize(n);
}

std::vector<node> SpSampler::randomPath() {
	node u = G.randomNode();
	node v = G.randomNode();
	while (u == v) {
		v = G.randomNode();
	}

	if (!G.isDirected() && cc.componentOfNode(u) != cc.componentOfNode(v)) {
		return std::vector<node>();
	}

	count endQ = 2;
	q[0] = u;
	q[1] = v;

	ballInd[u] = 1;
	ballInd[v] = 2;

	dist[u] = 0;
	dist[v] = 0;
	nPaths[u] = 1;
	nPaths[v] = 1;

	std::vector<std::pair<node, node>> spEdges;

	node x, randomEdge;
	bool hasToStop = false, useDegreeIn;
	count startU = 0, startV = 1, endU = 1, endV = 2, startCur, endCur,
	      *newEndCur;
	count sumDegsU = 0, sumDegsV = 0, *sumDegsCur;
	count totWeight = 0, curEdge = 0;

	auto procNeighbor = [&](node x, node y) {
		if (ballInd[y] == 0) {
			(*sumDegsCur) += getDegree(G, y, useDegreeIn);
			nPaths[y] = nPaths[x];
			ballInd[y] = ballInd[x];
			q[endQ++] = y;
			(*newEndCur)++;
			pred.addEdge(y, x);
			dist[y] = dist[x] + 1;
		} else if (ballInd[x] != ballInd[y]) {
			hasToStop = true;
			spEdges.push_back(std::make_pair(x, y));
		} else if (dist[y] == dist[x] + 1) {
			nPaths[y] += nPaths[x];
			pred.addEdge(y, x);
		}
	};

	while (!hasToStop) {
		if (sumDegsU <= sumDegsV) {
			startCur = startU;
			endCur = endU;
			startU = endQ;
			newEndCur = &endU;
			endU = endQ;
			sumDegsU = 0;
			sumDegsCur = &sumDegsU;
			useDegreeIn = false;
		} else {
			startCur = startV;
			endCur = endV;
			startV = endQ;
			newEndCur = &endV;
			endV = endQ;
			sumDegsV = 0;
			sumDegsCur = &sumDegsV;
			useDegreeIn = true;
		}

		while (startCur < endCur) {
			x = q[startCur++];

			if (useDegreeIn) {
				G.forInNeighborsOf(x, [&](node y) { procNeighbor(x, y); });
			} else {
				G.forNeighborsOf(x, [&](node y) { procNeighbor(x, y); });
			}
		}

		if (*sumDegsCur == 0) {
			hasToStop = true;
		}
	}

	if (spEdges.size() == 0) {
		removeAllEdges(endQ);
		std::fill(ballInd.begin(), ballInd.end(), 0);
		return std::vector<node>();
	}

	for (auto p : spEdges) {
		totWeight += nPaths[p.first] * nPaths[p.second];
	}

	randomEdge = Aux::Random::integer(totWeight - 1);
	std::vector<node> path;

	for (auto p : spEdges) {
		curEdge += nPaths[p.first] * nPaths[p.second];
		if (curEdge > randomEdge) {
			backtrackPath(u, v, p.first, path);
			backtrackPath(u, v, p.second, path);
			break;
		}
	}

	std::fill(ballInd.begin(), ballInd.end(), 0);
	removeAllEdges(endQ);
	return path;
}

void SpSampler::backtrackPath(const node u, const node v, const node start,
                              std::vector<node> &path) {
	if (start == u || start == v) {
		return;
	}

	count totWeight = nPaths[start];
	node randomPred, curPred = 0;
	node w = 0;

	path.push_back(start);
	randomPred = Aux::Random::integer(totWeight - 1);
	assert((pred.neighbors(start)).size() > 0);

	for (node t : pred.neighbors(start)) {
		w = t;
		curPred += nPaths[v];
		if (curPred > randomPred) {
			break;
		}
	}

	if (w != u && w != v) {
		backtrackPath(u, v, w, path);
	}
}

void SpSampler::removeAllEdges(const count endQ) {
	std::vector<node> resizedQ(endQ);
	std::copy(q.begin(), q.begin() + endQ, resizedQ.begin());
	pred.removeEdgesFromIsolatedSet(resizedQ);
}

count SpSampler::getDegree(const Graph &graph, node z, bool useDegreeIn) {
	return useDegreeIn ? graph.degreeIn(z) : graph.degree(z);
}
} // namespace NetworKit
