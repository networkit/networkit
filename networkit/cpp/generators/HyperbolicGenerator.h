/*
 * HyperbolicGenerator.h
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef HYPERBOLICGENERATOR_H_
#define HYPERBOLICGENERATOR_H_

#include <vector>
#include "../geometric/HyperbolicSpace.h"
#include "StaticGraphGenerator.h"
#include "../auxiliary/Timer.h"
#include "Quadtree/Quadtree.h"

namespace NetworKit {

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity
typedef index node; // node indices are 0-based

class HyperbolicGenerator: public NetworKit::StaticGraphGenerator {
public:

	/**
	 * @param[in] n Number of nodes
	 * @param[in] m Target number of edges
	 */
	HyperbolicGenerator(count n=10000, double avgDegree=6, double exp=3);


	Graph generate(const vector<double> &angles, const vector<double> &radii, Quadtree<index> &quad, double thresholdDistance);

	/**
	 * @param[in] angles Pointer to angles of node positions
	 * @param[in] radii Pointer to radii of node positions
	 * @param[in] r radius of poincare disk to place nodes in
	 * @param[in] thresholdDistance Edges are added for nodes closer to each other than this threshold
	 * @return Graph to be generated according to parameters
	 */
	Graph generate(const vector<double> &angles, const vector<double> &radii, double r, double thresholdDistance);
	Graph generateExternal(const vector<double> &angles, const vector<double> &radii, double k, double gamma);

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

	/**
	 * Set the capacity of a quadtree leaf.
	 *
	 * @param capacity Tuning parameter, default value is 1000
	 */
	void setLeafCapacity(count capacity) {
		this->capacity = capacity;
	}

	/**
	 * When using a theoretically optimal split, the quadtree will be flatter, but running time usually longer.
	 * @param split Whether to use the theoretically optimal split. Defaults to false
	 */
	void setTheoreticalSplit(bool split) {
		this->theoreticalSplit = split;
	}

	void setBalance(double balance) {
		this->balance = balance;
	}

	vector<double> getElapsedMilliseconds() {
		vector<double> result(threadtimers.size());
		for (index i = 0; i < result.size(); i++) {
			result[i] = threadtimers[i].elapsedMilliseconds();
		}
		return result;
	}

private:

	/**
	 * Set tuning parameters to their default values
	 */
	void initialize();

	/**
	 * graph parameters
	 */
	count nodeCount;
	double stretch;
	double factor;
	double alpha;

	/**
	 * tuning parameters
	 */
	count capacity;
	bool theoreticalSplit;
	double balance;
	static const bool directSwap = false;

	/**
	 * times
	 */
	vector<Aux::Timer> threadtimers;
};
}
#endif /* HYPERBOLICGENERATOR_H_ */
