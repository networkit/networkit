/*
 * AttributeGenerator.h
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#ifndef ATTRIBUTEGENERATOR_H_
#define ATTRIBUTEGENERATOR_H_

#include "../graph/Graph.h"
#include <unordered_map>
#include <vector>

namespace NetworKit {
/**
 * Abstract base class for graph attribute generator. It takes a graph (weighted or unweighted)
 * and calculates a graph attribute from the input graph.
 */
template<typename TInput, typename TOutput>
class AttributeGenerator {

public:

	/**
	 * Calculates an edge attribute for the edges of the given graph.
	 * (Possibly under consideration of the given attribute).
	 */
	virtual std::vector<TOutput> getAttribute(const Graph& g, const std::vector<TInput>& attribute) = 0;

	virtual ~AttributeGenerator() = default;

	std::vector<TOutput>* _getAttribute(Graph& g, std::vector<TInput>& attribute) {
		return new std::vector<TOutput>{std::move(getAttribute(g, attribute))};
	};

};

}


#endif /* ATTRIBUTEGENERATOR_H_ */
