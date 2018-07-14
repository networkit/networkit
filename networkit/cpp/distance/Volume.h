/*
 * Volume.h
 *
 * Created on: July 7, 2018
 * Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#ifndef VOLUME_H_
#define VOLUME_H_

#include "../graph/Graph.h"
#include <unordered_map>

namespace NetworKit {

/**
 * @ingroup volume
 */
class Volume {
private:
	/**
	 * For a given node n, the nodes within distance r are returned together with
	 * their distance to node n.
	 */
	static std::unordered_map<node, double> nodesWithinDistance(const Graph &G, double r, node n);

public:
	/**
	* Number of nodes within a given radius; average for many nodes
	*
	* @param G  the graph
	* @param r  the radius
	* @param samples	the number of samples
	*
	**/
	static double volume(const Graph &G, const double r, const count samples);

	/**
	* Number of nodes within different given radii; average for many nodes
	*
	* @param G  the graph
	* @param rs  the radii
	* @param samples	the number of samples
	*
	**/
	static std::vector<double> volume(const Graph &G, const std::vector<double> rs, const count samples);
};

} /* namespace NetworKit */

#endif /* VOLUME_H_ */
