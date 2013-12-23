/*
 * BetweennessCentrality.h
 *
 *  Created on: 09.12.2013
 *      Author: Lukas Barth, David Weiß
 */

#ifndef BetweennessCentrality_H_
#define BetweennessCentrality_H_

#include <list>

#ifdef __GNUC__
#include <ext/hash_map>
#else
#include <hash_map>
#endif

#include "../graph/Graph.h"

namespace std
{
 	using namespace __gnu_cxx;
}

namespace GrauBart {
using namespace NetworKit;

class BetweennessCentrality {

private:
	Graph *g;
	std::hash_map<node, std::vector<count>> dists; // This holds the BFS distance tables for every source vertex

	// This holds (for every source s) a map of vertices to the list of their *successors* (not predecessors).
	std::hash_map<node, std::hash_map<node, std::list<node>>> successors;

	// This holds (for every source s) the number of shortest paths from s to every vertex v
	std::hash_map<node, std::hash_map<node, count>> numShortestPaths;

	// This basically holds \delta_s for every vertex and is recursively built by computeDependencies()
	std::hash_map<node, std::hash_map<node, double>> dependencies;

	// This holds the final betweenness centralities
	std::hash_map<node, double> BCs;

	void computeTrees();
	void computeSuccessors();
	void computeDependencies();
	void accumulateDependencies();

public:
	BetweennessCentrality(Graph &G);

	virtual ~BetweennessCentrality();

	void run();

	std::hash_map<node, double> getCentrality();
};

} /* namespace NetworKit */
#endif /* BetweennessCentrality_H_ */
