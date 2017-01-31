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

	/**
	 * @param[in] bottomLeft Bottom left coordinate of the layout
	 * @param[in] topRight Top right coordinate of the layout
	 * @param[in] useGivenLayout Whether the current coordinates from the graph should be used
	 */
	Layouter(Point<float> bottomLeft, Point<float> topRight, bool useGivenLayout = false);

	virtual ~Layouter();

	/**
	 * Generates coordinates for the given graph g so that the graph can be drawn subsequently. 
	 * @param[in,out] g Graph to generate coordinates for
	 */
	virtual void draw(Graph& g) = 0;

	/**
	 * Assigns initial values to the coordinates.
	 * This is based on the given graph in case the given layout
	 * should be used. Otherwise the coordinates will be randomly initialized.
	 * @param[in] g Graph to initialize coordinates for
	 */
	virtual void initialize(Graph& g);

	/**
	 * Assigns random coordinates to every node in the given graph g within the given range on the layout.
	 * @param[in] g Graph to assign the coordinates to
	 */
	virtual void randomInitCoordinates(Graph& g);
};

} /* namespace NetworKit */
#endif /* LAYOUTER_H_ */
