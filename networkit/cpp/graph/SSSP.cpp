/*
 * SSSP.cpp
 *
 *  Created on: 15.04.2014
 *      Author: cls
 */

#include "SSSP.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

SSSP::SSSP(const Graph& G, node s, bool storePaths, bool storeStack, node target) : Algorithm(), G(G), source(s), target(target), storePaths(storePaths), storeStack(storeStack) {
}

std::vector<edgeweight> SSSP::getDistances(bool moveOut) {
	return (moveOut)?std::move(distances):distances;
}


std::vector<node> SSSP::getPath(node t, bool forward) const {
	if (! storePaths) {
		throw std::runtime_error("paths have not been stored");
	}
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


std::set<std::vector<node>> SSSP::getPaths(node t, bool forward) const {

	std::set<std::vector<node>> paths;
	if (previous[t].empty()) { // t is not reachable from source
		WARN("there is no path from ", source, " to ", t);
		return paths;
	}

	std::function<void (std::vector<node> suffix, node v) > trace = [&](std::vector<node> suffix, node v) {
		// base case
		suffix.push_back(v);
		if (v == source) {
			paths.insert(suffix);
			return;
		}
		for (node u : previous[v]) {
			trace(suffix, u);
		}
	};

	std::vector<node> emptyPath;
	trace(emptyPath, t);

	if (forward) {
		std::set<std::vector<node>> paths1;
		for (std::vector<node> path : paths) {
			std::reverse(std::begin(path), std::end(path));
			paths1.insert(path);
		}
		return paths1;
	}

	return paths;
}


std::vector<node> SSSP::getStack(bool moveOut) {
	if (!storeStack) {
		throw std::runtime_error("stack has not been stored");
	}
	return (moveOut)?std::move(stack):stack;
}

} /* namespace NetworKit */
