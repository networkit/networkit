/*
 * Layouter.cpp
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#include "../auxiliary/Random.h"

#include "Layouter.h"

namespace NetworKit {

Layouter::Layouter(Point<float> bottom_left, Point<float> top_right) :
		bottomLeft(bottom_left), topRight(top_right)
{
	// TODO Auto-generated constructor stub

}

Layouter::~Layouter() {
	// TODO Auto-generated destructor stub
}

Layouter::Layouter() {
}

void Layouter::randomInitCoordinates(Graph& g) {

	float x1 = bottomLeft.getValue(0);
	float y1 = bottomLeft.getValue(1);
	float width = topRight.getValue(0) - x1;
	float height = topRight.getValue(1) - y1;

	g.forNodes([&](node u) {
		float x = Aux::Random::real() * width + x1;
		float y = Aux::Random::real() * height + y1;

		TRACE("x: " , x);
		TRACE("y: " , y);

		std::vector<float> coords = {x, y};
		Point<float> p(coords);
		layout.push_back(p);
	});
}


} /* namespace NetworKit */
