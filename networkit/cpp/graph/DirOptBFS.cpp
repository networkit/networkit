/*
 * DirOptBFS.cpp
 *
 *  Created on: Feb 21, 2015
 *      Author: Maximilian Vogel
 */
#include <omp.h>
#include <queue>
#include <sstream>
#include "DirOptBFS.h"

namespace NetworKit {

DirOptBFS::DirOptBFS(const Graph& G, node source, bool storePaths, bool storeStack, count alpha, count beta, count max_threads) :
	SSSP(G, source, storePaths, storeStack, none),
	topdown(true),
	alpha(alpha),
	beta(beta),
	m_f(G.degree(source)),
	m_u((G.isDirected()?1:2) * G.numberOfEdges()),
	n_f(1),
	rhs_C_BT(G.numberOfNodes()/beta),
	max_threads((max_threads==0)?omp_get_max_threads():max_threads)
{
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges need to be indexed - call G.indexEdges() first");
	}
}

void DirOptBFS::run() {
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	count z = G.upperNodeIdBound();
	count currentDistance = 0;
	bool wasTopDown = true;
	count lastFrontierSize = 0;

	// prepare distance data structure
	distances.clear();
	distances.resize(z, infDist);

	// prepare data structures to store the paths
	if (storePaths) {
		previous.clear();
		previous.resize(z);
		npaths.clear();
		npaths.resize(z,0);
		previous[source] = {source};
		npaths[source] = 1;
	}

	// prepare the data structure to store the stack
	if (storeStack) {
		stack.reserve(G.numberOfNodes());
	}

	// prepare data structures to run the algorithm such as the visited flags for edges and nodes,
	// counters and queues
	frontier.clear();
	frontier.resize(z, 0);
	next.clear();
	next.resize(z, 0);
	std::vector<char> edgeActive(G.upperEdgeIdBound(),1);
	std::vector<char> visited(z, 0);
	std::vector<count> threadLocalCounter(max_threads, 0);

	// initialize the values for the source node
	qNext.push_back(source);
	distances[source] = currentDistance;
	visited[source] = 1;

	// stopping criterion: no new nodes have been added to the queue; depends on the step
	auto isFinished = [&]() {
		if (topdown) {
			return qNext.empty();
		} else {
			return n_f == 0;
		}
	};

	// implementation of section V of the paper.
	// note:
	//  - m_f is computed on the fly during the top-down step
	//  - m_u is computed during top-down step and, as it is easier, during prepareDatastructure instead of bottom-up step
	//  - n_f is not needed in the top-down step as it just the number of elements in the vector
	//  - n_f in bottom-up step is computed and reset in the else-case in the main loop.
	auto determineNextStep = [&]() {
		bool growing;
		if (topdown) {
			growing = qNext.size() >= lastFrontierSize;
			lastFrontierSize = qNext.size();
			topdown = !growing || m_f < (m_u / alpha);
		} else {
			growing = n_f >= lastFrontierSize;
			lastFrontierSize = n_f;
			topdown = !growing && n_f < rhs_C_BT;
		}
	};

	// implementation of the bottom-up step with paths 
	auto bottomUpStepWithPaths = [&](){
		G.balancedParallelForNodes([&](node v){
			G.forInNeighborsOf(v,[&](node, node u, edgeid eid){
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
			});
		});
	};

	// implementation of the bottom-up step
	auto bottomUpStep = [&](){
		G.balancedParallelForNodes([&](node v){
			if (!visited[v]) {
				bool found = false;
				G.forInNeighborsOf(v,[&](node, node u, edgeid eid) {
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
				});
			}
		});
	};


	// implementation of the top-down step
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

	// depending on the last step and the next step, data structures have to be swapped or converted
	auto prepareDatastructure = [&]() {
		if (topdown) {
			if (wasTopDown) {
				// qNext -> qFrontier
				std::swap(qFrontier,qNext);
			} else {
				// next -> qFrontier
				qFrontier.reserve(n_f);
				std::vector<std::vector<node>> threadLocalQueue(max_threads);
				#pragma omp parallel for schedule(guided)
				for (index i = 0; i < z; ++i) {
					if (next[i]) threadLocalQueue[omp_get_thread_num()].push_back(i);
				}
				for (const auto& q : threadLocalQueue) {
					qFrontier.insert(qFrontier.end(),q.begin(),q.end());
				}
				threadLocalCounter.assign(max_threads,0);
			}
			if (storeStack) {
				stack.insert(stack.end(),qFrontier.begin(),qFrontier.end());
			}
		} else {
			frontier.assign(z,false);
			if (wasTopDown) {
				// qNext -> frontier
				std::vector<count> mu(max_threads,0);
				#pragma omp parallel for schedule(guided)
				for (index i = 0; i < qNext.size(); ++i) {
					node current = qNext[i];
					frontier[current] = true;
					mu[omp_get_thread_num()] += G.degree(current);
				}
				for (const auto& m : mu) {
					m_u -= m;
				}
				if (storeStack) {
					stack.insert(stack.end(),qNext.begin(),qNext.end());
				}
				qNext.clear();
			} else {
				// next -> frontier
				std::swap(frontier, next);
				if (storeStack) {
					std::vector<std::vector<node>> threadLocalStack(max_threads);
					#pragma omp parallel for schedule(guided)
					for (index i = 0; i < z; ++i) {
						if (frontier[i]) {
							threadLocalStack[omp_get_thread_num()].push_back(i);
						}
					}
					for (const auto& s : threadLocalStack) {
						stack.insert(stack.end(),s.begin(),s.end());
					}
				}
			}
		}
	};

	// main loop
	while (!isFinished()) {
		++currentDistance;
		wasTopDown = topdown;
		determineNextStep();
		prepareDatastructure();
		if (topdown) {
			m_f = 0;
			topDownStep();
		} else {
			n_f = 0;
			if (storePaths)	bottomUpStepWithPaths();
			else bottomUpStep();
			count tmp = 0;
			#pragma omp parallel for reduction(+:tmp)
			for (index i = 0; i < max_threads; ++i) {
				tmp += threadLocalCounter[i];
			}
			threadLocalCounter.assign(max_threads,0);
			n_f = std::move(tmp);
		}
	}

	hasRun = true;
}

std::string DirOptBFS::toString() const {
	std::stringstream ss;
	ss << "DirOptBFS(storePaths=" << storePaths << ", storeStack=" << storeStack << ")";
	return ss.str();
}

} /* namespace NetworKit */
