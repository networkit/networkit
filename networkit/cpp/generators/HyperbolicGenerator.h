/*
 * HyperbolicGenerator.h
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef HYPERBOLICGENERATOR_H_
#define HYPERBOLICGENERATOR_H_

#include <vector>
#include <map>
#include "../geometric/HyperbolicSpace.h"
#include "StaticGraphGenerator.h"

namespace NetworKit {

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity
typedef index node; // node indices are 0-based

class HyperbolicGenerator: public NetworKit::StaticGraphGenerator {
public:
	HyperbolicGenerator();

	/**
	 * @param[in] n Number of nodes
	 * @param[in] factor Size of neighborhood radius. If factor=1, radius = R
	 * @param[in] alpha Dispersion parameter, default=1
	 * @param[in] stretchradius Stretching the hyperbolic disk results in thinner graphs, default=1
	 */
	HyperbolicGenerator(count n, double factor = 1, double alpha = 1, double stretchradius = 1);

	/**
	 * @param[in] n Number of nodes
	 * @param[in] m Number of edges
	 */
	HyperbolicGenerator(count n, count m);

	virtual ~HyperbolicGenerator();

	/**
	 * @param[in] n Number of nodes
	 * @param[in] stretch Parameter s
	 * @return Theoretical number of edges expected from a graph with alpha=1 and factor=1
	 */
	static double expectedNumberOfEdges(count n, double stretch);

	/**
	 * @param[in] angles Pointer to angles of node positions
	 * @param[in] radii Pointer to radii of node positions
	 * @param[in] R radius of hyperbolic disk to place nodes in
	 * @param[in] thresholdDistance Edges are added for nodes closer to each other than this threshold
	 * @return Graph to be generated according to parameters
	 */
	static Graph generate(vector<double> * angles, vector<double> * radii, double R, double thresholdDistance);

	/**
	 * Convenience function to convert polar coordinates into Cartesian coordinates
	 */
	static std::map<index, Point<float> > getCoordinates(vector<double> &angles, vector<double> &radii);

	/**
	 * @param[in] n Number of nodes
	 * @param[in] factor Size of neighborhood radius. If factor=1, radius = R
	 * @param[in] alpha Dispersion parameter, default=1
	 * @param[in] stretchradius Stretching the hyperbolic disk results in thinner graphs, default=1
	 * @return Graph to be generated according to parameters
	 */
	Graph generate(count n, double distanceFactor=1, double alpha=1, double stretchradius = 1);
	/**
	 * @return Graph to be generated according to parameters specified in constructor.
	 */
	Graph generate();

private:
	count nodeCount;
	double stretch;
	double factor;
	double alpha;

};
}
#endif /* HYPERBOLICGENERATOR_H_ */
