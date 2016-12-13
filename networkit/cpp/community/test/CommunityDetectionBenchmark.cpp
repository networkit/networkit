/*
 * CommunityDetectionBenchmark.h
 *
 *  Created on: 16.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#include <map>
#include <functional>

#include "CommunityDetectionBenchmark.h"
#include "../PLP.h"
#include "../PLM.h"
#include "../Modularity.h"
#include "../../centrality/Betweenness.h"
#include "../../centrality/PageRank.h"
#include "../../auxiliary/Timer.h"
#include "../../structures/Partition.h"

namespace NetworKit {

constexpr int runs = 20;

void CommunityDetectionBenchmark::SetUp() {

}

TEST_F(CommunityDetectionBenchmark, timeClusteringAlgos) {
	Aux::Timer timer;
	Modularity mod;

	// std::string graph = "../graphs/in-2004.graph";
	// std::string graph = "../graphs/uk-2002.graph";
	// std::string graph = "../graphs/uk-2007-05.graph";
	std::string graph = "input/polblogs.graph";

	printf("Reading graph file %s ...\n", graph.c_str());
	timer.start();
	const Graph G = this->metisReader.read(graph);
	timer.stop();
	printf("Reading graph took %.1f s\n", timer.elapsedMilliseconds() / 1000.0);

	for (int r = 0; r < runs; r++) {
		Graph Gcopy = G;
		PLP algo(Gcopy);

		timer.start();
		algo.run();
		Partition zeta = algo.getPartition();
		timer.stop();

		auto communitySizes = zeta.subsetSizes();


		printf("%s on %s: %.1f s\n\t# communities: %lu\n\tmodularity: %f\n",
			"Parallel Label Propagation", graph.c_str(),
			timer.elapsedMilliseconds() / 1000.0,
			zeta.numberOfSubsets(),
			mod.getQuality(zeta, G));
	}

	for (int r = 0; r < runs; r++) {
		Graph Gcopy = G;
		PLM algo(Gcopy);

		timer.start();
		algo.run();
		Partition zeta = algo.getPartition();
		timer.stop();

		auto communitySizes = zeta.subsetSizes();


		printf("%s on %s: %.1f s\n\t# communities: %lu\n\tmodularity: %f\n",
			"Parallel Louvain", graph.c_str(),
			timer.elapsedMilliseconds() / 1000.0,
			zeta.numberOfSubsets(),
			mod.getQuality(zeta, G));
	}
}

TEST_F(CommunityDetectionBenchmark, timePageRankCentrality) {
	Aux::Timer timer;

	// std::string graph = "../graphs/uk-2002.graph";
	std::string graph = "input/polblogs.graph";

const Graph G = this->metisReader.read(graph);

	for (int r = 0; r < runs; r++) {
		PageRank cen(G, 1e-6);

		timer.start();
		cen.run();
		timer.stop();
		auto ranking = cen.ranking();


		printf("%s on %s: %.1f s\n\tranking: [(%lu: %f), (%lu: %f), ...]\n",
			"Page Rank Centrality", graph.c_str(),
			timer.elapsedMilliseconds() / 1000.0,
			ranking[0].first, ranking[0].second,
			ranking[1].first, ranking[1].second);
	}
}

TEST_F(CommunityDetectionBenchmark, timeBetweennessCentrality) {
	Aux::Timer timer;

	// std::string graph = "../graphs/cond-mat-2005.graph";
	std::string graph = "input/polblogs.graph";

const Graph G = this->metisReader.read(graph);

	for (int r = 0; r < runs; r++) {
		Betweenness cen(G);

		timer.start();
		cen.run();
		timer.stop();
		auto ranking = cen.ranking();


		printf("%s on %s: %.1f s\n\tranking: [(%lu: %f), (%lu: %f), ...]\n",
			"Betweenness Centrality", graph.c_str(),
			timer.elapsedMilliseconds() / 1000.0,
			ranking[0].first, ranking[0].second,
			ranking[1].first, ranking[1].second);
	}
}


} /* namespace NetworKit */

#endif /*NOGTEST */
