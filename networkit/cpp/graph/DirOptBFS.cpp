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

DirOptBFS::DirOptBFS(const Graph& G, node source, bool storePaths, bool storeStack, count alpha, count beta, node target) :
	SSSP(G, source, storePaths, storeStack, target),
	hasQueuedNodes(false),
	topdown(true),
	alpha(alpha),
	beta(beta),
	m_f(0),
	m_u(0),
	n_f(0),
	rhs_C_BT(G.numberOfNodes()/beta) {
}

void DirOptBFS::run() {
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	count z = G.upperNodeIdBound();
	count currentDistance = 0;

	std::vector<bool> visited(z,false);

	auto wasVisited = [&](node v) {
		// numerous criteria could be used
		//return !previous[v].empty();
		//return distances[v] != infDist;
		return visited[v];
	};
	auto doTopDownStep = [&]() {
		//topdown = topdown?(m_f < (m_u/alpha)):(n_f < rhs_C_BT);
		if (topdown) {
			count m_f = 0;
			// manual computation of m_f if no bookkeeping is done
			// set would be better as one doesn't have to do a whole scan
			#pragma omp parallel for reduction(+:m_f)
			for (count u = 0; u < qFrontier.size(); ++u) {
				m_f += G.degree(u);
			}
			// manual computation of m_u if no bookkeeping is done
			count m_u = 0;
			#pragma omp parallel for reduction(+:m_u) 
			for (count u = 0; u < z; ++u) {
				m_u += (G.hasNode(u)&&!wasVisited(u))?G.degree(u):0;
			}
			//INFO("m_f > (m_u / alpha): ",m_f," / (",m_u," / ",alpha,")");
			topdown = !(m_f > (m_u / alpha));
		} else {
			//count n = G.numberOfNodes();
			count n_f = qFrontier.size();
			// manual computation of n_f is no bookkeeping is done
			/*#pragma omp parallel for reduction(+:n_f)
			for (count u = 0; u < qFrontier.size(); ++u) {
				n_f += frontier[u];
			}*/
			topdown = n_f < rhs_C_BT;
		}
		// reset bookkeeping variables that are computed for each iteration
		//m_f = 0;
		//n_f = 0;
		//INFO("topdown? ", topdown);
		return topdown;
	};


	auto relax = [&](node v) {
		if (!wasVisited(v)) {
			visited[v] = true;
			distances[v] = currentDistance;
			#pragma omp critical
			qNext.push_back(v);
		}
	};

	auto bottomUpStep = [&](){
		// this probably could be parallelized, however the following stuff needs to be sorted out
		// - concurrent writes on different indices vector<bool> and vector in general ? - Not a problem on general vectors (indices are different). But vector<bool>
		//    has a different implementation from all other vectors to save space, so this needs to be tested.
		// - synchronisation of bookkeeping variables m_f, n_f, m_u. The order of writes do not matter at all,
		//   however they need to be taken care of when parallelized. - These are the problem. Simply parallelizing this will result in
		//   wrong results due to race conditions. However, can the
		//G.forNodes([&](node v){
		G.balancedParallelForNodes([&](node v){
			// iterate over all nodes v that haven't been visited yet
			if (!wasVisited(v)) {
				// iterate over their neighbours u ...
				// G.forNeighborsOf is not an option because you can't "break" it
				for (auto &u : G.neighbors(v)) {
					// ... and if one of them belongs to the frontier ...
					if (frontier[u]) {
						// ... visit v
						relax(v);
						previous[v].push_back(u);
						// if we only want one predecessor/bfs tree, we break the loop here.
						//break;
					}
				}
			}
		});
	};

	auto topDownStep = [&]() {
		for (auto& current : qFrontier) {
			G.forNeighborsOf(current,[&](node v){
				relax(v);
				previous[v].push_back(current);
			});
		}
	};

	auto convertQueue = [&]() {
		for (auto& current : qFrontier) {
			frontier[current] = true;
		}
	};

	distances.clear();
	distances.resize(z, infDist);
	previous.clear();
	previous.resize(z);
	frontier.clear();
	frontier.resize(z,false);

	qFrontier.push_back(source);
	//hasQueuedNodes = true;
	previous[source] = {source};
	distances[source] = currentDistance;
	visited[source] = true;
/*	count time_topstep = 0;
	count time_botstep = 0;
	count n_top = 0;
	count n_bot = 0;
	Aux::Timer timer;*/
	//m_u = G.numberOfEdges() * 2;
	while (!qFrontier.empty()) {
	//while (hasQueuedNodes) {
		//hasQueuedNodes = false;
		//INFO("level: ",currentDistance, "\tsize of frontier: ",qFrontier.size());
		++currentDistance;
		if (doTopDownStep()) {
//			timer.start();
			topDownStep();
/*			timer.stop();
			time_topstep += timer.elapsedMilliseconds();
			n_top += n_f;*/
		} else {
			convertQueue();
//			timer.start();
			bottomUpStep();
/*			timer.stop();
			time_botstep += timer.elapsedMilliseconds();
			n_bot += n_f;*/
			// during a bottom-up step, we don't really handle the nodes in the frontier and clean it up
			// therefore, clear the frontier now. TODO: can this be done any cheaper?
			frontier.assign(z,false);
		}
		qFrontier.erase(qFrontier.begin(),qFrontier.end());
		std::swap(qFrontier,qNext);
	}
/*	INFO("time spent in top-down  steps:\t",time_topstep,"ms");
	INFO("nodes claimed in top-down  steps:\t",n_top);
	INFO("time spent in bottom-up steps:\t",time_botstep,"ms");
	INFO("nodes claimed in bottom-up steps:\t",n_bot);*/
}

} /* namespace NetworKit */
