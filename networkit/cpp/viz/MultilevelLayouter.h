/*
 * MultilevelLayouter.h
 *
 *  Created on: 27.01.2014
 *      Author: Henning
 */

#ifndef MULTILEVELLAYOUTER_H_
#define MULTILEVELLAYOUTER_H_

#include "Layouter.h"

namespace NetworKit {

/**
 * @ingroup viz
 */
// TODO: refactor to inherit from LayoutAlgorithm base class
class MultilevelLayouter: public NetworKit::Layouter {
protected:
	static const count N_THRSH;

public:
	MultilevelLayouter(Point<float> bottomLeft, Point<float> topRight, bool useGivenLayout = false);

	virtual void draw(Graph& G);
	virtual void drawInternal(Graph& G, count level);

	virtual void prolongCoordinates(Graph& Gcon, Graph& G);
};

} /* namespace NetworKit */
#endif /* MULTILEVELLAYOUTER_H_ */
