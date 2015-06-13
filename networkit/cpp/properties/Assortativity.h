/*
 * Assortativity.h
 *
 *  Created on: Jun 13, 2015
 *      Author: Christian Staudt
 */

#ifndef ASSORTATIVITY_H_
#define ASSORTATIVITY_H_

#include "../graph/Graph.h"


namespace NetworKit {

/**
 * @ingroup properties
 *
 * Assortativity computes a coefficient that expresses the correlation of a
 * node attribute among connected pairs of nodes.
 */
class Assortativity {

public:

	/**
	 * Create CoreDecomposition class for graph @a G.
	 *
	 * @param G The graph.
	 */
	Assortativity(const Graph& G, const std::vector<double>& attribute);

	/**
	*
	*/
	void run();


	/**
	* Return the assortativity coefficient.
	*/
	double getCoefficient() const;



private:

	const Graph& G;
	const std::vector<double>& attribute;
	bool ran; // whether algorithm has been run
};

} /* namespace NetworKit */
#endif /* ASSORTATIVITY_H_ */
