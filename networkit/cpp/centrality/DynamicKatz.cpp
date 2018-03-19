/*
 * DynamicKatz.cpp
 *
 *  Created on: 09.01.2015
 *      Author: Henning
 */

#include "DynamicKatz.h"
#include "../auxiliary/NumericTools.h"
#include <math.h>

namespace NetworKit {

DynamicKatz::DynamicKatz(const Graph& G, count k, bool groupOnly):
		Centrality(G, true), k(k), groupOnly(groupOnly)
{
	maxdeg = 0;
	G.forNodes([&](node u){
		if (G.degree(u) > maxdeg) {
			maxdeg = G.degree(u);
		}
	});
	alpha = double(1)/(maxdeg + 1);
	INFO("alpha: ", alpha);
	INFO("1/(1-alpha): ", 1/(1-alpha));
}

bool pairCompare(const std::pair<node, double>& firstElem, const std::pair<node, double>& secondElem) {
  return firstElem.second > secondElem.second;
}

void DynamicKatz::run() {
	INFO("Nodes: ", G.numberOfNodes(), ", edges: ", G.numberOfEdges());
	count z = G.upperNodeIdBound();
	count maxLevels = 15;
	nPaths.clear();
	nPaths.resize(maxLevels);
	for(count i = 0; i < maxLevels; i ++) {
		nPaths[i].resize(z, 1);
	}
//	std::vector<count> oldnPaths(z, 1);
	std::vector<bool> active(z, true);
	std::vector<double> ubound(z, 0.0);
	std::vector<double> lbound(z, 0.0);
	count nActive = z;

	scoreData.clear();
	scoreData.resize(z, 0.0);

	count totalActive = 0;
	count i = 1;
	count maxIter = 20;
	std::vector<std::pair<node, double>> orderedNodes;
	bool allOrdered = false; // required in the second phase. it is true when the first k nodes are ordered as well
	// while(nActive > k || !allOrdered) {
	while(nActive > k || (!groupOnly && !allOrdered)) {
		orderedNodes.clear();
		G.forNodes([&](node u){
				nPaths[i][u] = 0;
				G.forInEdgesOf(u, [&](node v, edgeweight ew) {
					nPaths[i][u] += nPaths[i-1][v];
				});
				scoreData[u] += pow(alpha, i) * nPaths[i][u];
				ubound[u] = scoreData[u] + nPaths[i][u] * pow(alpha, i + 1) * (1/(1 - alpha * maxdeg));
				lbound[u] = scoreData[u] + nPaths[i][u] * pow(alpha, i + 1) * (1/(1 - alpha));
				if(active[u]) {
					orderedNodes.push_back(std::make_pair(u, lbound[u]));
				}
		});
		// sort the nodes according to their lower bound
		std::sort(orderedNodes.begin(), orderedNodes.end(), pairCompare);


		INFO("Node with highest centrality: ", orderedNodes[0], ", score: ", scoreData[orderedNodes[0].first], ", ubound = ", ubound[orderedNodes[0].first]);
		INFO("K-th Node with highest centrality: ", orderedNodes[k-1], ", k = ", k);

		count removed = 0;
		for (count j = k; j < nActive; j ++) {
			node u = orderedNodes[j].first;
			assert(active[u]);
			//INFO("Node ", u, ", upper bound: ", ubound[u]);
			// if the upper bound on a node is smaller than the k-th lower bound, we can discard it
			if (ubound[u] <= orderedNodes[k - 1].second) {
				removed ++;
				active[u] = false;
				//INFO("Discarding ", u);
			}
		}
		nActive -= removed; //important to keep it, or you will not finish the loop!

		totalActive += nActive;
		// TODO careful: if two nodes have the same set of neighbors, they'll also have the same katz centrality!
		if (nActive == k) {
			allOrdered = true;
			// now we continue to iterate until the first k nodes are ordered
			count notOrdered = 0;
			for (count j = 1; j < k; j ++) {
				node first = orderedNodes[j-1].first, next = orderedNodes[j].first;
				if (ubound[next] > lbound[first]) {
					allOrdered = false;
					notOrdered ++;
				}
			}
			INFO("Not ordered: ", notOrdered);
		}
		i ++;
	}
	levelReached = i-1;
	INFO("Level reached: ", levelReached);

	// for a comparison, we run it until machine precision
	// nPaths.clear();
	// oldnPaths.clear();
	// active.clear();
	// ubound.clear();
	// lbound.clear();
	// nPaths.resize(z, 1);
	// oldnPaths.resize(z, 1);
	// active.resize(z, true);
	// ubound.resize(z, 0.0);
	// lbound.resize(z, 0.0);
	//
	// scoreData.clear();
	// scoreData.resize(z, 0.0);
	//
	// count changed = z;
	// i = 0;
	// while(changed > 0) {
	// 	changed = 0;
	// 	G.forNodes([&](node u){
	// 		// if (active[u]) { // TODO improve: iterate only over the active ones!
	// 			nPaths[u] = 0;
	// 			// note: inconsistency in definition in Newman's book (Ch. 7) regarding directed graphs
	// 			G.forInEdgesOf(u, [&](node v, edgeweight ew) {
	// 				nPaths[u] += oldnPaths[v];
	// 			});
	// 			edgeweight new_score = scoreData[u] + pow(alpha, i) * nPaths[u];
	// 			if (new_score > edgeweight(scoreData[u])) {
	// 				scoreData[u] = new_score;
	// 				changed ++;
	// 			}
	// 	});
	// 	G.forNodes([&](node u){
	// 			oldnPaths[u] = nPaths[u];
	// 	});
	// 	INFO("i = ", i, ", changed: ", changed);
	// 	i ++;
	// }
	// INFO("Iterations before reaching machine precision: ", i);

	hasRun = true;
}

void DynamicKatz::update(GraphEvent e){
	if (e.type != GraphEvent::EDGE_ADDITION && e.type != GraphEvent::EDGE_REMOVAL) {
		throw std::runtime_error("event type not allowed. Edge insertions or deletions only.");
	}
	if (e.type == GraphEvent::EDGE_REMOVAL && (nPaths[1][e.u] == 0 || nPaths[1][e.v] == 0)) {
		throw std::runtime_error("error: deleting an edge that did not exist before");
	}
	node u = e.u, v = e.v;
	count visitedEdges = 0;
	// first, we increase the Katz score of the two endpoints
	count z = G.upperNodeIdBound();
	std::vector<count> newPaths(z, 0);
	std::vector<count> newPathsPrevIt(z, 0);
	if (e.type == GraphEvent::EDGE_ADDITION) {
		scoreData[u] += alpha;
		scoreData[v] += alpha;
		newPathsPrevIt[u] = nPaths[1][u] + 1;
		newPathsPrevIt[v] = nPaths[1][v] + 1;
	} else {
		assert(nPaths[1][u] > 0);
		assert(nPaths[1][v] > 0);
		scoreData[u] -= alpha;
		scoreData[v] -= alpha;
		newPathsPrevIt[u] = nPaths[1][u] - 1;
		newPathsPrevIt[v] = nPaths[1][v] - 1;
	}
	std::vector<bool> isActive(z, false);
	std::vector<node> active;
	active.push_back(u);
	active.push_back(v);
	isActive[u] = true;
	isActive[v] = true;
	count i = 2;
	INFO("Starting update iteration");
	while(i <= levelReached) {
		// INFO("Level ", i, ", computing new paths");
		for (node u: active){
			newPaths[u] = nPaths[i][u];
			// INFO("u = ", u, ", npaths[", i, "] = ", nPaths[i][u]);
		}
		// INFO("Done initializing");
		std::queue<node> activated;
		for (node u: active){
			// notice: for directed graphs here the direction has to be the opposite of the static case
			G.forEdgesOf(u, [&](node v, edgeweight ew) {
				// INFO("Node ", u, ", neighbor ", v);
				visitedEdges ++;
				if(!isActive[v]) {
					isActive[v] = true;
					activated.push(v);
					newPaths[v] = nPaths[i][v];
				}
				// the old contrib should be subtracted only for the edges that existed before the insertion
				if ((v != e.v || u != e.u) && (v != e.u || u != e.v)) {
					newPaths[v] -= nPaths[i-1][u]; // subtract the old contribution and add the new one
				}
				newPaths[v] += newPathsPrevIt[u];
			});
			if (e.type == GraphEvent::EDGE_REMOVAL) {
				if (u == e.u) {
					newPaths[e.v] -= nPaths[i-1][u];
				} else if (u == e.v) {
					newPaths[e.u] -= nPaths[i-1][u];
				}
			}
		}
		INFO("Level ", i, ", updating scores");
		// now we update the scores (TODO and the bounds?)
		// also, we push the newly activated nodes into active
		while(activated.size() > 0) {
			node u = activated.front();
			activated.pop();
			active.push_back(u);
		}
		for (node u: active) {
			scoreData[u] -= pow(alpha, i) * nPaths[i][u];
			scoreData[u] += pow(alpha, i) * newPaths[u];
			if (newPathsPrevIt[u] > 0 || u == e.u || u == e.v) {
				nPaths[i-1][u] = newPathsPrevIt[u];
			} // otherwise I should not update this (this is a newly-discovered node)
			// Only for deletions, it might happen that u or v become isolated
			newPathsPrevIt[u] = newPaths[u];
		}
		INFO("i = ", i, ", nActive = ", active.size());
		i++;
	}
	INFO("Done update iteration. visitedEdges = ", visitedEdges, ", speedup: ", double(levelReached*G.numberOfEdges()*2)/visitedEdges);
	// last level: we need to update nPaths (TODO what about the bounds?)
	for (node u: active) {
		if (newPathsPrevIt[u] > 0 || u == e.u || u == e.v) {
			nPaths[i-1][u] = newPathsPrevIt[u];
		}
	}
	i --;
	// we compute the new bounds and check whether we still have the top-k nodes or we need more iterations
	std::vector<double> ubound(z);
	std::vector<double> lbound(z);
	std::vector<std::pair<node, double>> orderedNodes;
	G.forNodes([&](node u){
		ubound[u] = scoreData[u] + nPaths[i][u] * pow(alpha, i + 1) * (1/(1 - alpha * maxdeg));
		lbound[u] = scoreData[u] + nPaths[i][u] * pow(alpha, i + 1) * (1/(1 - alpha));
		orderedNodes.push_back(std::make_pair(u, lbound[u]));
	});
	std::sort(orderedNodes.begin(), orderedNodes.end(), pairCompare);
	bool allOrdered = true;
	for (count j = k; j < z; j ++) {
		node u = orderedNodes[j].first;
		if (ubound[u] > orderedNodes[k - 1].second) {
			allOrdered = false;
		}
	}
	INFO("All ordered? ", allOrdered);
	// if they are not all ordered, we might need new iterations... TODO
	// TODO what if the maxdeg increases???
}

} /* namespace NetworKit */
