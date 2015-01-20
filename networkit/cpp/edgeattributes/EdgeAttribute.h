/*
 * EdgeAttribute.h
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#ifndef EDGEATTRIBUTE_H_
#define EDGEATTRIBUTE_H_

#include "../graph/Graph.h"
#include <unordered_map>
#include <vector>

namespace NetworKit {
/**
 * Abstract base class for graph attribute generator.
 */
template<typename T>
class EdgeAttribute {

public:

	/**
	 * Calculates an edge attribute for the edges of the given graph.
	 */
	virtual std::vector<T> getAttribute() = 0;

	virtual ~EdgeAttribute() = default;

	std::vector<T>* _getAttribute() {
		return new std::vector<T>{std::move(getAttribute())};
	};

};

}


#endif /* EDGEATTRIBUTE_H_ */
