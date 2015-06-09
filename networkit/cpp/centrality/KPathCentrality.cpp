/*
 * KPathCentrality.cpp
 *
 *  Created on: 05.10.2014
 *      Author: nemes
 */

#include <stack>

#include "KPathCentrality.h"
#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

KPathCentrality::KPathCentrality(const Graph& G, double alpha, count k) : Centrality(G, false, false) {
	if (alpha >= -0.5 && alpha <= 0.5) {
		this->alpha = alpha;
	} else {
		throw std::runtime_error("alpha must lie in interval [-0.5, 0.5]");
	}
	if (k == 0) {
		this->k = log(G.numberOfNodes() + G.numberOfEdges());
	} else if (k >0) {
		this->k = k;
	} else {
		throw std::runtime_error("k must be an integer");
	}
}

void KPathCentrality::run() {
	count z = G.upperNodeIdBound();
	count n = G.numberOfNodes();
	scoreData.clear();
	scoreData.resize(z);

	std::vector<count> counter;
	std::vector<bool> explored;

	counter.assign(z, 0);
	explored.assign(z, false);

	count t = 2 * k * k * pow(n, 1 - 2 * alpha) * log(n);
	std::stack<node> stack;
	node v;

	for (index i = 1; i <= t; i++) { // FIXME: int -> count
		node s = G.randomNode();
		auto l = Aux::Random::integer(1, k);
		explored[s] = true;
		stack.push(s);
		count j = 1;

		while (j <= l) {
			edgeweight sum = 0;
			std::vector<node> neighbours;
			neighbours.clear();
			std::vector<edgeweight> weights;
			weights.clear();
			G.forNeighborsOf(s, [&](node u, edgeweight ew) {
				if (!explored[u]) {
					neighbours.push_back(u);
					weights.push_back(1/ew);
					sum += 1/ew;
				}
			});
			if (neighbours.empty()) {
				break;
			}
			if (G.isWeighted()) {
				double random = Aux::Random::real(0, sum);
				for (index x = 0; x < weights.size(); x++) {
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

	G.forNodes([&](node v) {
		scoreData[v] = k * n * ((double) counter[v] / t);
	});

	hasRun = true;
}


} /* namespace NetworKit */
