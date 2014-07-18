/*
 * SSSP.cpp
 *
 *  Created on: 15.04.2014
 *      Author: cls
 */

#include "SSSP.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

SSSP::SSSP(const Graph& G, node s) : G(G), source(s) {
	
}

std::vector<edgeweight> SSSP::getDistances() const {
	return distances;
}


std::vector<node> SSSP::getPath(node t, bool forward) const {
	std::vector<node> path;
	if (previous[t].empty()) { // t is not reachable from source
		WARN("there is no path from ", source, " to ", t);
		return path;
	}
	node v = t;
	while (v != source) {
		path.push_back(v);
		v = previous[v].front();
	}
	path.push_back(source);

	if (forward) {
		std::reverse(path.begin(), path.end());
	}
	return path;
}


std::set<std::vector<node> > SSSP::getPaths(node t, bool forward) const {
	throw std::runtime_error("FIXME: correct implementation needed");

	std::set<std::vector<node> > paths;
	if (previous[t].empty()) { // t is not reachable from source
		WARN("there is no path from ", source, " to ", t);
		return paths;
	}


	std::function<std::set<std::vector<node> > (std::vector<node>& prefix, node v) > trace = [&](std::vector<node>& prefix, node v) {
		// base case

		prefix.push_back(v);
		std::set<std::vector<node> > paths;
		paths.insert(prefix);
		for (node u : previous[v]) {
			auto returned = trace(prefix, u);
			paths.insert(returned.begin(), returned.end());
		}
		return paths;
	};

	std::vector<node> emptyPath;
	auto thePaths = trace(emptyPath, t);

	if (forward) {
		for (auto path : thePaths) {
			std::reverse(path.begin(), path.end());
		}
	}
	return thePaths;
}

} /* namespace NetworKit */
