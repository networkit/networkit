
#include "clique.h"

namespace NetworKit {

Clique::Clique(const Graph& G, const unsigned int missingEdges) : G(G), missingEdges(missingEdges) {
}

std::vector<std::set<node> > Clique::run(node& seed) {
	std::set<std::set<node> > oldsets;
	std::set<std::set<node> > newsets;

	auto unset_intersect = [&] (const std::set<node>& seta, const std::set<node>& setb) {
		std::set<node> result;
		for (auto a : seta) {
			if (setb.find(a) != setb.end()) {
				result.insert(a);
			}
		}
		return result;
	};

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

		for (const auto& oldset : oldsets) {
			for (auto u : oldset) {
				for (auto v : G.neighbors(u)) {
					if (oldset.find(v) != oldset.end())
						break;
					auto nv = G.neighbors(v);
					auto intersection = unset_intersect(oldset, std::set<node>(nv.begin(), nv.end()));
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
	if (missingEdges > 0) {

		for (const auto& oldset : oldsets) {
			for (const auto u : oldset) {
				for (auto v : G.neighbors(u)) {
					if (oldset.find(v) != oldset.end())
						break;
					auto nv = G.neighbors(v);
					auto intersection = unset_intersect(oldset, std::set<node>(nv.begin(), nv.end()));
					if (intersection.size() >= oldset.size() - missingEdges) {
						auto baum = oldset;
						baum.insert(v);
						newsets.insert(baum);
					}
				}
			}
		}
		if (!newsets.empty()) {
			oldsets = newsets;
		}
	}

	std::vector<std::set<node> > result(oldsets.begin(), oldsets.end());
	return result;
}

}
