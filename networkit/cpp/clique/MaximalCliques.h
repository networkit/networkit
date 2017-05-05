#ifndef CLIQUE_H_
#define CLIQUE_H_

#include "../graph/Graph.h"
#include "../base/Algorithm.h"

#include <unordered_map>

namespace NetworKit {

class MaximalCliques : public Algorithm {

public:
	MaximalCliques(const Graph& G);

	void run() override;

	const std::vector<std::vector<node>>& getCliques() const;

protected:
	const Graph& G;

	std::vector<std::vector<node>> result;

	struct StaticOutGraph;

	void tomita(const StaticOutGraph& outGraph, std::vector<node>& pxvector, std::vector<index>& pxlookup, index xbound, index xpbound, index pbound, std::vector<node>& r);

	node findPivot(const StaticOutGraph& outGraph, std::vector<node>& pxvector, std::vector<index>& pxlookup, index xbound, index xpbound, index pbound);

	struct StaticOutGraph {
		StaticOutGraph(const Graph& G, const std::vector<index>& pxlookup) : firstOut(G.upperNodeIdBound() + 1), head(G.numberOfEdges()) {
			index currentOut = 0;
			for (node u = 0; u < G.upperNodeIdBound(); ++u) {
				firstOut[u] = currentOut;
				if (G.hasNode(u)) {
					index xpboundU = pxlookup[u];
					G.forEdgesOf(u, [&](node v) {
						if (xpboundU < pxlookup[v]) {
							head[currentOut++] = v;
						}
					});
				}
			}
			firstOut[G.upperNodeIdBound()] = currentOut;

		}

		template <typename F>
		void forOutEdgesOf(node u, F callback) const {
			for (index i = firstOut[u]; i < firstOut[u + 1]; ++i) {
				callback(head[i]);
			}
		}

		bool hasNeighbor(node u, node v) const {
			for (index i = firstOut[u]; i < firstOut[u + 1]; ++i) {
				if (head[i] == v) return true;
			}

			return false;
		}

		std::vector<index> firstOut;
		std::vector<node> head;
	};
};

}

#endif /* CLIQUE_H_ */
