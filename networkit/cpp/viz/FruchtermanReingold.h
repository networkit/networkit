/*
 * FruchtermanReingold.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#ifndef FRUCHTERMANREINGOLD_H_
#define FRUCHTERMANREINGOLD_H_

#include "Layouter.h"
#include "Point.h"

#include <vector>
#include <cmath>


namespace NetworKit {

const count MAX_ITER = 300;
const double EPS = 0.1;

/**
 * @ingroup viz
 * Fruchterman-Reingold graph drawing algorithm. We mostly follow
 * the description in Stephen G. Kobourov: Spring Embedders and Force
 * Directed Graph Drawing Algorithms.
 */

// TODO: refactor to inherit from LayoutAlgorithm base class
class FruchtermanReingold: public NetworKit::Layouter {
private:
	static const float INITIAL_STEP_LENGTH;
	static const float OPT_PAIR_SQR_DIST_SCALE;

	count maxIter;
	float prec;
	float step;

public:

	/**
	 * Constructor. DO NOT use for creating objects, only needed for
	 * Python shell.
	 */
	FruchtermanReingold() {}

	/**
	 * Constructor.
	 * @param[in] bottomLeft Coordinate of point in bottom/left corner
	 * @param[in] topRight Coordinate of point in top/right corner
	 */
	FruchtermanReingold(Point<float> bottomLeft, Point<float> topRight, bool useGivenCoordinates = false, count maxIterations = MAX_ITER, float precision = EPS);

	/**
	 * Assigns coordinates to vertices in graph @a g
	 * @param[in,out] g Graph used to provide the vertices
	 */
	virtual void draw(Graph& g);
};

} /* namespace NetworKit */
#endif /* FRUCHTERMANREINGOLD_H_ */
