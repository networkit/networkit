/*
 * ChibaNishizekiTriangleCounter.h
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#ifndef CHIBANISHIZEKI_H_
#define CHIBANISHIZEKI_H_

#include "../graph/Graph.h"
#include "AttributeGenerator.h"

namespace NetworKit {

/**
 * An implementation of the triangle counting algorithm by Chiba/Nishizeki.
 */
class ChibaNishizekiTriangleCounter : public AttributeGenerator<int, int> {

public:

	std::vector<int> getAttribute(const Graph& graph, const std::vector<int>& attribute);
	~ChibaNishizekiTriangleCounter() = default;

private:
	void removeNode(Graph& graph, node u);
};

} /* namespace NetworKit */

#endif /* CHIBANISHIZEKI_H_ */
