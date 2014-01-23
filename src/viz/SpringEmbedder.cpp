/*
 * SpringEmbedder.cpp
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#include "SpringEmbedder.h"

namespace NetworKit {

SpringEmbedder::SpringEmbedder(Point<float> bottom_left, Point<float> top_right):
		Layouter(bottom_left, top_right)
{

}

SpringEmbedder::SpringEmbedder() {
}

SpringEmbedder::~SpringEmbedder() {

}

} /* namespace NetworKit */
