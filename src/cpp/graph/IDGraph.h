/*
 * IDGraph.h
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef IDGRAPH_H
#define IDGRAPH_H

#include <limits>
#include <string>
#include <algorithm>

#include "IGraph.h"
#include "../Globals.h"
#include "../viz/Point.h"

namespace NetworKit {

/**
 * Interface for all directed graph classes. Every graph class has to implement all interface methods.
 */
class IDGraph : public virtual IGraph {

public:

	struct NodeDegree {
		count in;
		count out;
		NodeDegree():
			in(0),
			out(0)
		{}
		count total() const { return in + out; }
	};

	// TODO: structure for Edge (start, end, weight)?


	/** NODE PROPERTIES **/

	/**
	 * Return the number of neighbors for node v.
	 */
	virtual NodeDegree degree(node v) const = 0;

	/**
	 * Return the number of incoming edges to node v. For an undirected graph this is the
	 * same as degree(v).
	 */
	count degreeIn(node v) const { return this->degree(v).in; }

	/**
	 * Return the number of outgoing edges from node v. For an undirected graph this is the
	 * same as degree(v).
	 */
	count degreeOut(node v) const { return this->degree(v).out; }
	
	/**
	 * @return In-Volume of the node, which is the sum of all incoming edges.
	 */
	virtual edgeweight volume(node v) const;
};





} /* namespace NetworKit */

#endif /* IDGRAPH_H */
