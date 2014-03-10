/*
 * Betweenness.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "Betweenness.h"
#include "../auxiliary/PrioQueue.h"
#include <stack>
#include <queue>

namespace NetworKit {

Betweenness::Betweenness(const Graph& G, bool normalized) : Centrality(G, normalized) {

}

void Betweenness::run() {
	using namespace std;
	count z = G.upperNodeIdBound();
	scoreData.clear();
	scoreData.resize(z);

	// TODO: reduce code duplication; not entirely avoidable due to different data types
	if (G.isWeighted()) {
		static const edgeweight INF = numeric_limits<edgeweight>::max();
		G.forNodes([&] (node s) {
			/* Nodes in order of increasing distance from s. */
			stack<node> increasing;

			/* Determine the shortest path DAG from s via Dijkstra's algorithms. */
			vector<set<pair<edgeweight, node> > > parents(z); /* For each vertex pair of distance, node of parents in shortest-path DAG. */
			vector<count> nshort(z);            /* Number of shortest paths. */
			vector<edgeweight> lshort(z, INF);  /* Length of shortest paths (= INF means unvisited). */

			// Dijkstra's algorithm that also computes the shortest path DAG
			lshort[s] = 0;
			nshort[s] = 1;
			Aux::PrioQueue<edgeweight, node> pq(lshort); // priority queue with distance-node pairs

			auto relax([&](node u, node v, edgeweight ew)
			{
				if (lshort[v] > lshort[u] + ew) { // part of shortest path so far
					// record length of shortest path
					lshort[v] = lshort[u] + ew;

					// adapt PQ entry
					pq.decreaseKey(lshort[v], v);
				}

				if (lshort[v] >= lshort[u] + ew) {
					// shortest path of at least equal length => record
					nshort[v] += nshort[u];

					while (parents[v].size() > 0) {
						std::pair<edgeweight, node> distNodePair = (* parents[v].begin());

						if (distNodePair.first > lshort[v]) {
							// wrong parent distNodePair.second => subtract its number of shortest paths
							nshort[v] -= nshort[distNodePair.second];
							parents[v].erase(distNodePair);
						}
						else break;
					}

					parents[v].insert({lshort[v], u});
				}
			});

			while (pq.size() > 0) {
				node v = pq.extractMin().second;
				increasing.push(v);

				G.forWeightedNeighborsOf(v, [&] (node w, edgeweight ew) {
					relax(v, w, ew);
				});
			}

			/* Now compute the dependencies in order of decreasing distance. */
//			DEBUG("Processing DAG of node ", s);
			vector<double> dependency(z, 0);
			while (! increasing.empty()) {
				node w = increasing.top();
				increasing.pop();
//				DEBUG("DAG node ", w, "; number of shortest paths: ", nshort[w]);

				for (auto entry: parents[w]) {
					/* Recursive formula: see lecture. */
					node v = entry.second;
					dependency[v] += double(nshort[v])/nshort[w] * (1 + dependency[w]);
				}
				if (w != s) {
					scoreData[w] += dependency[w];
				}
			}
		});
	}
	else {
		static const count INF = numeric_limits<count>::max();
		G.forNodes([&] (node s) {
			/* Nodes in order of increasing distance from s. */
			stack<node> increasing;

			/* Determine the shortest path tree from s via BFS. */
			vector<vector<node>> parents(z);          /* Parents in DAG. */
			vector<count> nshort(z), lshort(z, INF);  /* Number and length of shortest paths (= INF means unvisited). */
			queue<node> bfs_queue;                    /* Working queue. */

			/* BFS that also computes the shortest path DAG. */
			bfs_queue.push(s);
			lshort[s] = 0;
			nshort[s] = 1;
			while (!bfs_queue.empty()) {
				node v = bfs_queue.front();
				bfs_queue.pop();
				increasing.push(v);

				G.forNeighborsOf(v, [&] (node w) {
					/* w found for first time -> enqueue */
					if (lshort[w] == INF) {
						bfs_queue.push(w);
						lshort[w] = lshort[v] + 1;
					}

					/* Another shortest path to w via v. */
					if (lshort[w] == lshort[v] + 1) {
						nshort[w] = nshort[w] + nshort[v];
						parents[w].push_back(v);
					}
				});
			}

			/* Now compute the dependencies in order of decreasing distance. */
			vector<double> dependency(z, 0);
			while (!increasing.empty()) {
				node w = increasing.top();
				increasing.pop();

				for (node v: parents[w]) {
					/* Recursive formula: see lecture. */
					dependency[v] += double(nshort[v])/nshort[w] * (1 + dependency[w]);
				}
				if (w != s) {
					scoreData[w] += dependency[w];
				}
			}
		});
	}

	if (normalized) {
		// divide by the number of possible pairs
		count n = G.numberOfNodes();
		count pairs = (n-2) * (n-1);
		G.parallelForNodes([&](node u){
			scoreData[u] = scoreData[u] / pairs;
		});
	}

}






} /* namespace NetworKit */

