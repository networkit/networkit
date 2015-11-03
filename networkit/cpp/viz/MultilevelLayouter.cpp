/*
 * MultilevelLayouter.cpp
 *
 *  Created on: 27.01.2014
 *      Author: Henning
 */

#include "MultilevelLayouter.h"
#include "FruchtermanReingold.h"
#include "MaxentStress.h"
#include "../community/PLP.h"
#include "../io/METISGraphWriter.h"
#include "../community/EdgeCut.h"
#include "../coarsening/ParallelPartitionCoarsening.h"

namespace NetworKit {

const count MultilevelLayouter::N_THRSH = 15;

MultilevelLayouter::MultilevelLayouter(Point<float> bottomLeft, Point<float> topRight, bool useGivenLayout):
		Layouter(bottomLeft, topRight, useGivenLayout)
{

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
		METISGraphWriter gWriter;
		gWriter.write(G, true, "output/test-multi-coarsest.graph");
		DEBUG("initial layout by FR, G's size: ", G.numberOfNodes());
		FruchtermanReingold initLayouter(bottomLeft, topRight, false);
		initLayouter.draw(G);
		PostscriptWriter writer;
		writer.write(G, "output/test-multi-coarsest-FR.eps");
		MaxentStress layouter(bottomLeft, topRight, true);
		layouter.draw(G);
		writer.write(G, "output/test-multi-coarsest-ME.eps");
	}
	else {
		// compute clustering
		PLP clusterer(G);
		clusterer.run();
		Partition clustering = clusterer.getPartition();
		EdgeCut ec;
		INFO("Clustering: #clusters: ", clustering.numberOfSubsets(), "; cut: ", ec.getQuality(clustering, G));

		// coarsen by clustering
		ParallelPartitionCoarsening contracter(G, clustering);
		contracter.run();
		Graph Gcon = contracter.getCoarseGraph();
		std::vector<node> mapping = contracter.getFineToCoarseNodeMapping();

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
		FruchtermanReingold layouter(bottomLeft, topRight, true); //, 50 * (level + 1), 0.1); // TODO: externalize
		layouter.draw(G);
	}
}

} /* namespace NetworKit */
