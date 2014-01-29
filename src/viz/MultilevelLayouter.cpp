/*
 * MultilevelLayouter.cpp
 *
 *  Created on: 27.01.2014
 *      Author: Henning
 */

#include "MultilevelLayouter.h"

namespace NetworKit {

const count MultilevelLayouter::N_THRSH = 15;

MultilevelLayouter::MultilevelLayouter(Point<float> bottomLeft, Point<float> topRight, bool useGivenLayout):
		Layouter(bottomLeft, topRight, useGivenLayout)
{

}

MultilevelLayouter::~MultilevelLayouter() {

}

void MultilevelLayouter::prolongCoordinates(Graph& Gcon, Graph& G) {

}

void MultilevelLayouter::draw(Graph& G) {
	drawWrapper(G, 0);
}

void MultilevelLayouter::drawWrapper(Graph& G, count level) {
	count n = G.numberOfNodes();

	if (n <= N_THRSH) {
		// unrecursive part: call drawing routine
		FruchtermanReingold layouter(bottomLeft, topRight, false);
		layouter.draw(G);
	}
	else {
		// compute matching
		LocalMaxMatcher matcher(none);
		Matching M = matcher.run(G);

		// coarsen by matching
		MatchingContracter contracter;
		auto mypair = contracter.run(G, M, true);
		Graph& Gcon = mypair.first;
		auto mapping = mypair.second;

		// make recursive call
		drawWrapper(Gcon, level + 1);

		// apply recursive solution to current graph
		G.initCoordinates();
		G.forNodes([&](node v) {
			G.setCoordinate(v, Gcon.getCoordinate(mapping[v]));
//			TRACE("coordinate of ", v, ": ", G.getCoordinate(v, 0), " / ", G.getCoordinate(v, 1));
		});
		DEBUG("local refinement of graph of size ", n);

		// run drawing code on current graph
		FruchtermanReingold layouter(bottomLeft, topRight, true, 150 * (level + 1), 0.1); // TODO: externalize
		layouter.draw(G);
	}
}

} /* namespace NetworKit */
