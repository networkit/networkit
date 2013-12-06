/*
 * Centrality_Hoske.cpp
 *
 *  Created on: 05.12.2013
 *      Author: dhoske
 */

#include "Centrality_Hoske.h"

namespace NetworKit {

/** Computes the betweenness centrality of the nodes in G. */
std::vector<double> betweennessCentrality_Hoske(const Graph& G) {
	using namespace std;
	static const count INF = numeric_limits<count>::max();
	count n = G.numberOfNodes();
	vector<double> betweenness(n);

	G.forNodes([&] (node s) {
		/* Nodes in order of increasing distance from s. */
		stack<node> increasing;

		/* Determine the shortest path dag from s via bfs. */
		vector<vector<node>> parents(n);         /* Parents in dag. */
		vector<count> nshort(n), lshort(n, INF);  /* Number and length of shortest paths. */
		queue<node> bfs_queue;                   /* Working queue. */

		bfs_queue.push(s);
		lshort[s] = 0; /* lshort < 0 marks a vertex that has not been visited! */
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
		vector<double> dependency(n);
		while (!increasing.empty()) {
			node w = increasing.top(); increasing.pop();
			//cout << "Node: " << w << " " << s << endl;
			for (node v: parents[w]) {
				/* Recursive formula: see lecture. */
				dependency[v] += double(nshort[v])/nshort[w] * (1 + dependency[w]);
			}
			if (w != s) {
				betweenness[w] += dependency[w];
			}
		}
	});

	return betweenness;
}

} /* namespace NetworKit */
