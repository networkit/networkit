/*
 * Betweenness.cpp
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#include "Betweenness.h"
#include <stack>
#include <queue>

namespace NetworKit {

Betweenness::Betweenness(const Graph& G) : Centrality(G) {

}

void Betweenness::run() {
	using namespace std;
	static const count INF = numeric_limits<count>::max();
	count z = G.upperNodeIdBound();
	scoreData.clear();
	scoreData.resize(z);

	G.forNodes([&] (node s) {
		/* Nodes in order of increasing distance from s. */
		stack<node> increasing;

		/* Determine the shortest path dag from s via bfs. */
		vector<vector<node>> parents(z);          /* Parents in dag. */
		vector<count> nshort(z), lshort(z, INF);  /* Number and length of shortest paths (= INF means unvisited). */
		queue<node> bfs_queue;                    /* Working queue. */

       /* BFS that also computes the shortest path dag. */
		bfs_queue.push(s);
		lshort[s] = 0;
		nshort[s] = 1;
		while (!bfs_queue.empty()) {
			node v = bfs_queue.front(); bfs_queue.pop();
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
		vector<double> dependency(z);
		while (!increasing.empty()) {
			node w = increasing.top(); increasing.pop();
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






} /* namespace NetworKit */

