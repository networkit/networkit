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

	edgeAttribute getAttribute(const Graph& graph);

private:
	void triangleFound(edgeAttribute& map, const node& u, const node& v, const node& w);
	void removeNode(Graph& graph, const node& u);
};

} /* namespace NetworKit */

#endif /* CHIBANISHIZEKI_H_ */
