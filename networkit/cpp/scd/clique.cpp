
#include "clique.h"

#include <algorithm>

namespace NetworKit {

Clique::Clique(const Graph& G, const unsigned int missingEdges) : G(G), missingEdges(missingEdges) {
}

std::vector<std::set<node> > Clique::run(node& seed) {
	std::set<std::set<node> > oldsets;
	std::set<std::set<node> > newsets;

	// create first set of (not maximum) cliques
	for (auto u : G.neighbors(seed)) {
		std::set<node> baum;
		baum.insert(seed);
		baum.insert(u);
		newsets.insert(baum);
	}

	// search for maximum cliques
	while (!newsets.empty()) {
		oldsets = newsets;
		newsets = std::set<std::set<node> >();

		for (auto oldset : oldsets) {
			for (auto u : oldset) {
				for (auto v : G.neighbors(u)) {
					if (oldset.find(v) != oldset.end())
						break;
					auto nu = G.neighbors(v);
					std::set<node> intersection;
					std::set_intersection(oldset.begin(), oldset.end(), nu.begin(), nu.end(), std::inserter(intersection, intersection.begin()));
					if (intersection.size() == oldset.size()) {
						auto baum = oldset;
						baum.insert(v);
						newsets.insert(baum);
					}
				}
			}
		}
	}

	// extend search allowing 'missingEdges'
	newsets = oldsets;
	while (!newsets.empty()) {
		oldsets = newsets;
		newsets = std::set<std::set<node> >();

		for (auto oldset : oldsets) {
			for (auto u : oldset) {
				for (auto v : G.neighbors(u)) {
					if (oldset.find(v) != oldset.end())
						break;
					auto nu = G.neighbors(v);
					std::set<node> intersection;
					std::set_intersection(oldset.begin(), oldset.end(), nu.begin(), nu.end(), std::inserter(intersection, intersection.begin()));
					if (intersection.size() >= oldset.size() - missingEdges) {
						auto baum = oldset;
						baum.insert(v);
						newsets.insert(baum);
					}
				}
			}
		}
	}

	std::vector<std::set<node> > result(oldsets.begin(), oldsets.end());
	return result;
}

}
