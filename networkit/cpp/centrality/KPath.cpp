/*
 * KPath.cpp
 *
 *  Created on: 05.10.2014
 *      Author: nemes
 */

#include <stack>

#include "KPath.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

KPath::KPath(const Graph& G, bool normalized) : Centrality(G, normalized, false) {

}

void KPath::run() {
	//TODO parameters
	double alpha = 0.3;
	count k = 5;
	count z = G.upperNodeIdBound();
	count n = G.numberOfNodes();
	scoreData.clear();
	scoreData.resize(z);

	std::vector<count> counter;
	std::vector<bool> explored;

	counter.assign(z, 0);
	explored.assign(z, false);

	count t = 2 * k * k * pow(n, 1 - 2 * alpha) * log2(n);
	std::stack<node> stack;
	node v;

	for (int i = 1; i <= t; i++) {
		node s = G.randomNode();
		int l = Aux::Random::integer(1, k);
		explored[s] = true;
		stack.push(s);
		int j = 1;
		std::vector<node> neighbours;

		while (j <= l) {
			neighbours.clear();
			G.forNeighborsOf(s, [&](node u) {
				if (!explored[u]) {
					neighbours.push_back(u);
				}
			});
			if (neighbours.empty()) {
				break;
			}
			if (G.isWeighted()) {
				// TODO
				edgeweight sum = 0;
				std::vector<edgeweight> weights;
				weights.clear();

				G.forNeighborsOf(s, [&](node u, edgeweight ew) {
					if (!explored[u]) {
						weights.push_back(1/ew);
						sum += 1/ew;
					}
				});
				double random = Aux::Random::real(0, sum);
				for (int x = 0; x < weights.size(); x++) {
					if (random < weights[x]) {
						v = neighbours[x];
						break;
					}
					random -= weights[x];
				}
			} else {
				v = neighbours[Aux::Random::integer(0,neighbours.size() - 1)];
			}
			explored[v] = true;
			stack.push(v);
			counter[v]++;
			s = v;
			j++;
		}
		while (!stack.empty()) {
			v = stack.top();
			stack.pop();
			explored[v] = false;
		}
	}

	G.forNodes(
			[&](node v) {
				std::cout << v << ": " << k << " " << n << " " << " " << counter[v] << " " << t << std::endl;
				scoreData[v] = k * n * ((double) counter[v] / t);
			});

	if (normalized) {
		// TODO
	}

}


} /* namespace NetworKit */
