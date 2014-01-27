/*
 * MultilevelLayouter.cpp
 *
 *  Created on: 27.01.2014
 *      Author: Henning
 */

#include "MultilevelLayouter.h"

namespace NetworKit {

const count MultilevelLayouter::N_THRSH = 15;

MultilevelLayouter::MultilevelLayouter(Point<float> bottomLeft, Point<float> topRight):
		Layouter(bottomLeft, topRight)
{

}

MultilevelLayouter::~MultilevelLayouter() {

}

void MultilevelLayouter::prolongCoordinates(Graph& Gcon, Graph& G) {

}

void MultilevelLayouter::draw(Graph& G) {
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
		std::pair<Graph, NodeMap<node> > mypair = contracter.run(G, M, true);
		Graph& Gcon = mypair.first;
		NodeMap<node>& mapping = mypair.second;

		// make recursive call
		draw(Gcon);

		// apply recursive solution to current graph
		G.initCoordinates();
		G.forNodes([&](node v) {
			G.setCoordinate(v, 0, Gcon.getCoordinate(mapping[v], 0));
			G.setCoordinate(v, 1, Gcon.getCoordinate(mapping[v], 1));

			DEBUG("coordinate of ", v, ": ", G.getCoordinate(v, 0), " / ", G.getCoordinate(v, 1));
		});
		DEBUG("local refinement of graph of size ", n);

		// run drawing code on current graph
		FruchtermanReingold layouter(bottomLeft, topRight, true);
		layouter.draw(G);
	}

}

} /* namespace NetworKit */
