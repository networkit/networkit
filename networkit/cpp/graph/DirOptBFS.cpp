/*
 * DirOptBFS.cpp
 *
 *  Created on: Feb 21, 2015
 *      Author: Maximilian Vogel
 */

#include <queue>
#include "DirOptBFS.h"
#include "../auxiliary/Timer.h"

namespace NetworKit {

DirOptBFS::DirOptBFS(const Graph& G, node source, bool storePaths, bool storeStack, count alpha, count beta) : 
	SSSP(G, source, storePaths, storeStack), 
	hasQueuedNodes(false), 
	topdown(true), 
	alpha(alpha), 
	beta(beta), 
	m_f(0), 
	m_u(0), 
	n_f(0), 
	rhs_C_BT(G.numberOfNodes()/beta) {
}

void DirOptBFS::run(node t) {
	count currentDistance = 0;
	auto doTopDownStep = [&]() {
		topdown = topdown?(m_f < (m_u/alpha)):(n_f < rhs_C_BT);
		/*if (topdown) {
			count m_f = 0;
			// manual computation of m_f is no bookkeeping is done
			// set would be better as one doesn't have to do a whole scan
			#pragma omp parallel for reduction(+:m_f)
			for (count u = 0; u < frontier.size(); ++u) {
				if (frontier[u]) m_f += G.degree(u);
			}
			// manual computation of m_u is no bookkeeping is done
			count m_u = G.parallelSumForNodes([&](node v){
				return (previous[v].empty()?G.degree(v):0);
			});
			topdown = !(m_f > (m_u / alpha));
		} else {
			count n = G.numberOfNodes();
			count n_f = 0;
			// manual computation of n_f is no bookkeeping is done
			#pragma omp parallel for reduction(+:n_f)
			for (count u = 0; u < frontier.size(); ++u) {
				n_f += frontier[u];
			}
			topdown = n_f < rhs_C_BT;
		}*/
		// reset bookkeeping variables that are computed for each iteration
		m_f = 0;
		n_f = 0;
		return topdown;
	};


	auto bottomUpStep = [&](){
		// this probably could be parallelized, however the following stuff needs to be sorted out
		// - concurrent writes on different indices vector<bool> and vector in general ?
		// - synchronisation of bookkeeping variables m_f, n_f, m_u. The order of writes do not matter at all, 
		//   however they need to be taken care of when parallelized.
		G.forNodes([&](node v){
			// iterate over all nodes v that haven't been visited yet
			if (previous[v].empty()) {
				// iterate over their neighbours u ...
				for (auto &u : G.neighbors(v)) {
					// ... and if one of them belongs to the frontier ...
					if (frontier[u]) {
						// ... visit v
						previous[v].push_back(u);
						distances[v] = currentDistance + 1;
						// add it to the frontier for the next iteration.
						next[v] = true;
						// do some "bookkeeping" for the direction optimizing heuristic
						m_f += G.degree(v);
						n_f++;
						hasQueuedNodes = true;
						m_u -= G.degree(u);
						break; // if we only want one predecessor/bfs tree, we break the loop here.
					}
				}
			}
		});
	};

	auto topDownStep = [&]() {
		// TODO:
		// this scan could only be avoided with a queue or set-like data structure
		// however, this will degrade performance as the current frontier[u]-accesses in the bottom-up-step
		// will result in "frontier.find(u) != frontier.end()" which seems to be way more expensive
		for (count u = 0, end = frontier.size(); u < end; ++u) {
			// if the node is not in the frontier, just continue
			if(!frontier[u]) continue;
			// unset the node in the frontier
			frontier[u] = false;
			G.forNeighborsOf(u,[&](node v) {
				// if the node hasn't been visited yet...
				if (previous[v].empty()) {
					// ... visit it
					previous[v].push_back(u);
					distances[v] = currentDistance + 1;
					// add it to the frontier for the next iteration.
					next[v] = true;
					// do some "bookkeeping" for the direction optimizing heuristic
					n_f++;
					m_f += G.degree(v);
					hasQueuedNodes = true;

				}
			});
			m_u -= G.degree(u);
		}
	};

	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	count z = G.upperNodeIdBound();
	distances.clear();
	distances.resize(z, infDist);
	previous.clear();
	previous.resize(z);
	frontier.clear();
	frontier.resize(z,false);
	next.clear();
	next.resize(z,false);

/*	std::vector<bool> visited;
	visited.resize(z, false);*/
	frontier[source] = true;
	hasQueuedNodes = true;
	previous[source] = {source};
	distances[source] = currentDistance;
/*	count time_topstep = 0;
	count time_botstep = 0;
	count n_top = 0;
	count n_bot = 0;
	Aux::Timer timer;*/
	m_u = G.numberOfEdges() * 2;
	while (hasQueuedNodes) {
		hasQueuedNodes = false;
		if (doTopDownStep()) {
//			timer.start();
			topDownStep();
/*			timer.stop();
			time_topstep += timer.elapsedMilliseconds();
			n_top += n_f;*/
		} else {
//			timer.start();
			bottomUpStep();
/*			timer.stop();
			time_botstep += timer.elapsedMilliseconds();
			n_bot += n_f;*/
			// during a bottom-up step, we don't really handle the nodes in the frontier and clean it up
			// therefore, clear the frontier now. TODO: can this be done any cheaper?
			frontier.assign(z,false);
		}
		++currentDistance;
		std::swap(frontier,next);
	}
/*	INFO("time spent in top-down  steps:\t",time_topstep,"ms");
	INFO("nodes claimed in top-down  steps:\t",n_top);
	INFO("time spent in bottom-up steps:\t",time_botstep,"ms");
	INFO("nodes claimed in bottom-up steps:\t",n_bot);*/
}

} /* namespace NetworKit */
