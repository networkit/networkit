/*
 * ForestFireAttributizer.cpp
 *
 *  Created on: 26.08.2014
 *      Author: Gerd Lindner
 */

#include "ForestFireAttributizer.h"
#include <limits>
#include <set>
#include <queue>
#include "../auxiliary/Log.h"

namespace NetworKit {

ForestFireAttributizer::ForestFireAttributizer(double pf, double targetBurntRatio): pf(pf), targetBurntRatio(targetBurntRatio) {}

std::vector<double> ForestFireAttributizer::getAttribute(const Graph& graph, const std::vector<int>& attribute) {

	std::vector<count> burnt (graph.upperEdgeIdBound(), 0);
	count edgesBurnt = 0;

	while (edgesBurnt < targetBurntRatio * graph.numberOfEdges()) {
		//Start a new fire
		std::queue<node> activeNodes;
		std::vector<bool> visited (graph.upperNodeIdBound(), false);
		activeNodes.push(graph.randomNode());

		auto forwardNeighbors = [&](node u) {
			std::vector<node> validEdges;
			graph.forNeighborsOf(u, [&](node x){
				if (! visited[x]) {
					validEdges.push_back(x);
				}
			});
			return validEdges;
		};

		while (! activeNodes.empty()) {
			node v = activeNodes.front();
			activeNodes.pop();

			std::vector<node> validNeighbors = forwardNeighbors(v);
			std::set<node> burntNeighbors;
			while (true) {
				double q = Aux::Random::real(1.0);
				if (q > pf || validNeighbors.empty()) {
					break;
				}
				count index = Aux::Random::integer(validNeighbors.size() - 1);
				burntNeighbors.insert(validNeighbors[index]);
				validNeighbors[index] = validNeighbors.back();
				validNeighbors.pop_back();
			}

			for (node x : burntNeighbors) {
				activeNodes.push(x);
				burnt[graph.edgeId(v, x)]++;
				edgesBurnt++;
				visited[x] = true;
			}
		}
	}

	std::vector<double> burntNormalized (graph.numberOfEdges(), 0.0);
	double maxv = (double) *std::max_element(std::begin(burnt), std::end(burnt));

	if (maxv > 0) {
		count idx = 0;
		for (auto& b : burnt) {
			burntNormalized[idx++] = b / maxv;
		}
	}

	return burntNormalized;
}

} /* namespace NetworKit */
