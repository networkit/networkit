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
class ChibaNishizekiTriangleCounter : public AttributeGenerator {

public:

	EdgeAttribute getAttribute(const Graph& graph, const EdgeAttribute& attribute);

private:
	void triangleFound(const Graph& graph, EdgeAttribute& map, node u, node v, node w);
	void removeNode(Graph& graph, node u);
};

} /* namespace NetworKit */

#endif /* CHIBANISHIZEKI_H_ */
