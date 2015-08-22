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

DirOptBFS::DirOptBFS(const Graph& G, node source, count alpha, count beta, node target) :
	SSSP(G, source, false, false, target),
	topdown(true),
	alpha(alpha),
	beta(beta),
	m_f(G.degree(source)),
	m_u(0),
	n_f(1),
	rhs_C_BT(G.numberOfNodes()/beta) {
}

void DirOptBFS::run() {
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	count z = G.upperNodeIdBound();
	count currentDistance = 0;

	distances.clear();
	distances.resize(z, infDist);
	//previous.clear();
	//previous.resize(z);
	frontier.clear();
	frontier.reserve(z);
	//frontier.resize(z,false);
	std::vector<bool> visited(z,false);
	count max_threads = omp_get_max_threads();
	std::vector<std::vector<node>> threadLocalNext(max_threads);

	qNext.push_back(source);
	//previous[source] = {source};
	distances[source] = currentDistance;
	visited[source] = true;
	bool wasTopDown = true;
	count lastFrontierSize = 0;
	bool growing;

	auto wasVisited = [&](node v) {
		// numerous criteria could be used
		//return !previous[v].empty();
		//return distances[v] != infDist;
		return visited[v];
	};

	auto isFinished = [&]() {
		if (topdown) {
			return qNext.empty();
		} else {
			for (auto &q : threadLocalNext) {
				if (!q.empty()) return false;
			}
			return true;
		}
	};

	auto determineNextStep = [&]() {
		if (topdown) {
			growing = qNext.size() > lastFrontierSize;
			lastFrontierSize = qNext.size();
			// manual computation of m_u
			// TODO: can this be computed on the fly?
			count m_u = 0;
			#pragma omp parallel for reduction(+:m_u) 
			for (count u = 0; u < z; ++u) {
				m_u += (G.hasNode(u)&&!wasVisited(u))?G.degree(u):0;
			}
			topdown = m_f < (m_u / alpha) && growing;
		} else {
			for (auto& q : threadLocalNext) {
				n_f += q.size();
			}
			growing = n_f > lastFrontierSize;
			lastFrontierSize = n_f;
			topdown = n_f < rhs_C_BT && !growing;
		}
	};

	auto bottomUpStep = [&](){
		G.balancedParallelForNodes([&](node v){
			// iterate over all nodes v that haven't been visited yet
			if (!wasVisited(v)) {
				// iterate over their neighbours u ...
				// G.forNeighborsOf is not an option because you can't "break" it
				for (auto &u : G.neighbors(v)) {
					// ... and if one of them belongs to the frontier ...
					if (frontier[u]) {
						// ... visit v
						if (!wasVisited(v)) {
							visited[v] = true;
							distances[v] = currentDistance;
							count tid = omp_get_thread_num();
							threadLocalNext[tid].push_back(v);
						}
						//previous[v].push_back(u);
						// if we only want one predecessor/bfs tree, we break the loop here.
						break;
					}
				}
			}
		});
	};

	auto relax = [&](node v) {
		if (!wasVisited(v)) {
			visited[v] = true;
			distances[v] = currentDistance;
			qNext.push_back(v);
			m_f += G.degree(v);
		}
	};

	auto topDownStep = [&]() {
		for (auto& current : qFrontier) {
			G.forNeighborsOf(current,[&](node v){
				relax(v);
				//previous[v].push_back(current);
			});
		}
	};

	auto prepareDatastructure = [&]() {
		if (topdown) {
			qFrontier.clear();
			if (wasTopDown) {
				// qNext -> qFrontier
				std::swap(qFrontier,qNext);
			} else {
				// threadLocalNext -> qFrontier
				qFrontier.reserve(n_f);
				for (auto& q : threadLocalNext) {
					qFrontier.insert(qFrontier.end(),q.begin(),q.end());
				}
				threadLocalNext.clear();
				threadLocalNext.resize(max_threads);
			}
		} else {
			frontier.assign(z,false);
			if (wasTopDown) {
				// qNext -> frontier
				for (auto& current : qNext) {
					frontier[current] = true;
				}
				qNext.clear();
			} else {
				// threadLocalNext -> frontier
				for (auto& q : threadLocalNext) {
					for (auto& current : q) {
						frontier[current] = true;
					}
				}
				threadLocalNext.clear();
				threadLocalNext.resize(max_threads);
			}
		}
	};

/*	count time_topstep = 0;
	count time_botstep = 0;
	count n_top = 0;
	count n_bot = 0;
	Aux::Timer timer;*/
	//m_u = G.numberOfEdges() * 2;
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
			bottomUpStep();
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
