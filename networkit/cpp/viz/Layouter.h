/*
 * Layouter.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#ifndef LAYOUTER_H_
#define LAYOUTER_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * DEPRECATED: use base class LayoutAlgorithm instead
 * @ingroup viz
 */
class Layouter {
protected:
	Point<float> bottomLeft;
	Point<float> topRight;
	std::vector<Point<float> > layout;
	bool initNecessary;

public:
	/**
	 * DO NOT use to construct objects. Nullary constructor needed for Python shell.
	 */
	Layouter() {}
	Layouter(Point<float> bottomLeft, Point<float> topRight, bool useGivenLayout = false);
	virtual ~Layouter();

	virtual void draw(Graph& g) = 0;

	virtual void initialize(Graph& g);
	virtual void randomInitCoordinates(Graph& g);
};

} /* namespace NetworKit */
#endif /* LAYOUTER_H_ */
