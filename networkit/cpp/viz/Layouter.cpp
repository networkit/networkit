/*
 * Layouter.cpp
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#include "../auxiliary/Random.h"

#include "Layouter.h"

namespace NetworKit {

Layouter::Layouter(Point<float> bottom_left, Point<float> top_right, bool useGivenLayout) :
		bottomLeft(bottom_left), topRight(top_right), initNecessary(! useGivenLayout)
{

}

Layouter::~Layouter() {

}

void Layouter::initialize(Graph& G) {
	if (initNecessary) {
		// initialize randomly
		randomInitCoordinates(G);
		G.initCoordinates();
	}
	else {
		// feed coordinates from g into layout
		layout.resize(G.numberOfNodes());
		G.parallelForNodes([&](node v) {
			layout[v] = G.getCoordinate(v);
		});
	}
}

void Layouter::randomInitCoordinates(Graph& g) {

	float x1 = bottomLeft[0];
	float y1 = bottomLeft[1];
	float width = topRight[0] - x1;
	float height = topRight[1] - y1;

	layout.resize(g.numberOfNodes());

	g.parallelForNodes([&](node u) {
		float x = Aux::Random::real() * width + x1;
		float y = Aux::Random::real() * height + y1;

//		TRACE("x: " , x);
//		TRACE("y: " , y);

		Point<float> p(x, y);
		layout[u] = p;
	});
}


} /* namespace NetworKit */
