/*
 * BackboneCalculator.h
 *
 *  Created on: 21.05.2014
 *      Author: Gerd Lindner
 */

#ifndef BACKBONECALCULATOR_H_
#define BACKBONECALCULATOR_H_

#include "../graph/Graph.h"
#include "AttributeGenerator.h"

namespace NetworKit {

/** 
 * Abstract base class for Backbone Calculators.
 */
class BackboneCalculator {

public:
	/**
	 * Calculates the backbone graph for the given input graph.
	 */
	virtual Graph calculate(const Graph& g, const edgeAttribute& attribute) = 0;

	/** only to be used by cython - this eliminates an unnecessary copy */
	Graph* _calculate(const Graph& g, const edgeAttribute& attribute) {
		return new Graph{std::move(calculate(g, attribute))};
	};

	/**
	 * Creates a new undirected graph that contains only the nodes of the given graph.
	 * TODO: Implement a clone method in Graph instead?
	 * @param graph 	the original graph to copy
	 * @param weighted	whether the new graph should be weighted
	 */
	Graph cloneNodes(const Graph& graph, bool weighted);
};

} /* namespace NetworKit */
#endif /* BACKBONECALCULATOR_H_ */
