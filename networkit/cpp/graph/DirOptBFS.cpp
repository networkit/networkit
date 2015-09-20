/*
 * DirOptBFS.cpp
 *
 *  Created on: Feb 21, 2015
 *      Author: Maximilian Vogel
 */
#include <omp.h>
#include <queue>
#include "DirOptBFS.h"
#include "../auxiliary/Timer.h"

namespace NetworKit {

DirOptBFS::DirOptBFS(const Graph& G, node source, bool storePaths, bool storeStack, count alpha, count beta, node target) :
	SSSP(G, source, storePaths, storeStack, target),
	topdown(true),
	alpha(alpha),
	beta(beta),
	m_f(G.degree(source)),
	m_u((G.isDirected()?1:2) * G.numberOfEdges()),
	n_f(1),
	rhs_C_BT(G.numberOfNodes()/beta) {
}

void DirOptBFS::run() {
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	count z = G.upperNodeIdBound();
	count currentDistance = 0;

	distances.clear();
	distances.resize(z, infDist);
	std::vector<char> edgeActive;
	if (G.hasEdgeIds())
		edgeActive.resize(G.upperEdgeIdBound(),1);
	if (storePaths) {
		previous.clear();
		previous.resize(z);
		npaths.clear();
		npaths.resize(z,0);
		previous[source] = {source};
		npaths[source] = 1;
	}

	if (storeStack) {
		std::stack<node> empty;
		std::swap(stack, empty);
	}
	frontier.clear();
	frontier.resize(z, 0);
	next.clear();
	next.resize(z, 0);
	std::vector<char> visited(z, 0);
	count max_threads = omp_get_max_threads();
	std::vector<count> threadLocalCounter(max_threads, 0);

	qNext.push_back(source);
	distances[source] = currentDistance;
	visited[source] = 1;

	bool wasTopDown = true;
	count lastFrontierSize = 0;
	bool growing;

	auto isFinished = [&]() {
		if (topdown) {
			return qNext.empty();
		} else {
			for (auto &s : threadLocalCounter) {
				if (s > 0) return false;
			}
			return true;
		}
	};

	auto determineNextStep = [&]() {
		if (topdown) {
			growing = qNext.size() >= lastFrontierSize;
			lastFrontierSize = qNext.size();
			topdown = !growing || m_f < (m_u / alpha);
		} else {
			count tmp = 0;
			#pragma omp parallel for reduction(+:tmp)
			for (index i = 0; i < max_threads; ++i) {
				tmp += threadLocalCounter[i];
			}
			n_f = std::move(tmp);
			growing = n_f >= lastFrontierSize;
			lastFrontierSize = n_f;
			topdown = !growing && n_f < rhs_C_BT;
		}
	};

	auto bottomUpStepWithPaths = [&](){
		G.balancedParallelForNodes([&](node v){
			G.forInNeighborsOf(v,[&](node v, node u, edgeid eid){
			//for (auto &u : G.inNeighbors(v)) {
				if (edgeActive[eid] && frontier[u]) {
					edgeActive[eid] = 0;
					if (!visited[v]) {
						visited[v] = 1;
						distances[v] = currentDistance;
						count tid = omp_get_thread_num();
						threadLocalCounter[tid] += 1;
						next[v] = 1;
						previous[v] = {u};
						npaths[v] = npaths[u];
					} else if (distances[v]-1 == distances[u]) {
						previous[v].push_back(u);
						npaths[v] += npaths[u];
					}
				}
			//}
			});
		});
	};

	auto bottomUpStep = [&](){
		G.balancedParallelForNodes([&](node v){
			if (!visited[v]) {
				bool found = false;
				G.forInNeighborsOf(v,[&](node v, node u, edgeid eid) {
				//for (auto &u : G.inNeighbors(v)) {
					if (edgeActive[eid] && frontier[u]) {
						edgeActive[eid] = 0;
						if (!found) {
							visited[v] = 1;
							distances[v] = currentDistance;
							count tid = omp_get_thread_num();
							threadLocalCounter[tid] += 1;
							next[v] = 1;
							found = true;
						}
					}
				//}
				});
			}
		});
	};


	auto topDownStep = [&]() {
		while (!qFrontier.empty()) {
			auto current = qFrontier.back();
			qFrontier.pop_back();

			G.forNeighborsOf(current,[&](node, node v, edgeid eid){
				if (edgeActive[eid]) {
					edgeActive[eid] = 0;
					if (!visited[v]) {
						visited[v] = 1;
						distances[v] = currentDistance;
						qNext.push_back(v);
						m_f += G.degree(v);
						if (storePaths) {
							previous[v] = {current};
							npaths[v] = npaths[current];
						}
					} else if (storePaths && distances[v]-1 == distances[current]) {
							previous[v].push_back(current);
							npaths[v] += npaths[current];
					}
				}
			});
			m_u -= G.degree(current);
		}
	};

	auto prepareDatastructure = [&]() {
		if (topdown) {
			if (wasTopDown) {
				// qNext -> qFrontier
				std::swap(qFrontier,qNext);
			} else {
				// next -> qFrontier
				for (index i = 0; i < z; ++i) {
					if (next[i]) qFrontier.push_back(i);
				}
				threadLocalCounter.assign(max_threads,0);
			}
			if (storeStack) {
				for (auto& u : qFrontier) {
					stack.push(u);
				}
			}
		} else {
			frontier.assign(z,false);
			if (wasTopDown) {
				// qNext -> frontier
				for (auto& current : qNext) {
					frontier[current] = true;
					m_u -= G.degree(current);
					if (storeStack) {
							stack.push(current);
					}
				}
				qNext.clear();
			} else {
				// next -> frontier
				if (storeStack) {
					for (index i = 0; i < z; ++i) {
						if (next[i]) {
							stack.push(i);
							frontier[i] = 1;
						} else {
							frontier[i] = 0;
						}
						next[i] = 0;
					}
				} else {
					std::swap(frontier, next);
					next.assign(z,0);
				}
				threadLocalCounter.assign(max_threads,0);
			}
		}
	};

/*	count time_topstep = 0;
	count time_botstep = 0;
	count n_top = 0;
	count n_bot = 0;
	Aux::Timer timer;*/
	while (!isFinished()) {
		++currentDistance;
		wasTopDown = topdown;
		determineNextStep();
		prepareDatastructure();
		if (topdown) {
			m_f = 0;
//			timer.start();
			topDownStep();
/*			timer.stop();
			time_topstep += timer.elapsedMilliseconds();
			n_top += n_f;*/
		} else {
			n_f = 0;
//			timer.start();
			if (storePaths)	bottomUpStepWithPaths();
			else bottomUpStep();
/*			timer.stop();
			time_botstep += timer.elapsedMilliseconds();
			n_bot += n_f;*/
		}
	}
/*	INFO("time spent in top-down  steps:\t",time_topstep,"ms");
	INFO("nodes claimed in top-down  steps:\t",n_top);
	INFO("time spent in bottom-up steps:\t",time_botstep,"ms");
	INFO("nodes claimed in bottom-up steps:\t",n_bot);*/
}

} /* namespace NetworKit */
