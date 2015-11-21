/*
 * Assortativity.h
 *
 *  Created on: Jun 13, 2015
 *      Author: Christian Staudt
 */

#ifndef ASSORTATIVITY_H_
#define ASSORTATIVITY_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"
#include "../structures/Partition.h"


namespace NetworKit {

/**
 * @ingroup correlation
 *
 * Assortativity computes a coefficient that expresses the correlation of a
 * node attribute among connected pairs of nodes.
 */
class Assortativity : public Algorithm {

public:

	/**
	 * Initialize Assortativity with a graph @a G and an array of numerical
	 * node values.
	 *
	 * @param G The graph.
	 * @param attribute		numerical node value array
	 */
	Assortativity(const Graph& G, const std::vector<double>& attribute);


	/**
	* Initialize Assortativity with a graph @a G and a partition of the node set
	*
	* @param G The graph.
	* @param partition		partition of the node set
	*/
	Assortativity(const Graph& G, const Partition& partition);

	/**
	*
	*/
	void run() override;


	/**
	* Return the assortativity coefficient.
	*/
	double getCoefficient() const;

	/**
	 * Returns a string with the algorithm's name and its parameters, if there are any. Subclasses should override it.
	 * @return The string representation of the algorithm.
	 */
	std::string toString() const override;

	bool isParallel() const override;


private:
	const Graph& G;
	const std::vector<double> emptyVector;
	const Partition emptyPartition;
	const std::vector<double>& attribute;
	const Partition& partition;
	bool nominal; // whether we calculate assortativity for a nominal or ordinal attribute
	double coefficient;
};

} /* namespace NetworKit */
#endif /* ASSORTATIVITY_H_ */
