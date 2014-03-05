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
	drawInternal(G, 0);
}

void MultilevelLayouter::drawInternal(Graph& G, count level) {
	count n = G.numberOfNodes();

	if (n <= N_THRSH) {
		// unrecursive part: call drawing routine
		FruchtermanReingold layouter(bottomLeft, topRight, false);
		layouter.draw(G);
		PostscriptWriter writer(G);
		writer.write("output/test-multi-coarsest.eps");
	}
	else {
		// compute clustering
		PLP clusterer;
		Partition clustering = clusterer.run(G);

		// coarsen by clustering
		ClusterContracter contracter;
		auto mypair = contracter.run(G, clustering);
		Graph& Gcon = mypair.first;
		auto mapping = mypair.second;

		// make recursive call
		drawInternal(Gcon, level + 1);

		// apply recursive solution to current graph
		G.initCoordinates();
		G.forNodes([&](node v) {
			G.setCoordinate(v, Gcon.getCoordinate(mapping[v]));
//			TRACE("coordinate of ", v, ": ", G.getCoordinate(v, 0), " / ", G.getCoordinate(v, 1));
		});
		DEBUG("local refinement of graph of size ", n);

		// run drawing code on current graph
		FruchtermanReingold layouter(bottomLeft, topRight, true, 50 * (level + 1), 0.1); // TODO: externalize
		layouter.draw(G);
	}
}

} /* namespace NetworKit */
