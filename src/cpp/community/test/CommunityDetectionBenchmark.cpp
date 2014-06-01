/*
 * CommunityDetectionBenchmark.h
 *
 *  Created on: 16.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#include <map>

#include "CommunityDetectionBenchmark.h"
#include "../PLP.h"
#include "../PLM.h"
#include "../Modularity.h"
#include "../../centrality/Betweenness.h"
#include "../../centrality/PageRank.h"
#include "../../auxiliary/Timer.h"
#include "../../structures/Partition.h"

namespace NetworKit {

void CommunityDetectionBenchmark::SetUp() {

}

TEST_F(CommunityDetectionBenchmark, timeClusteringAlgos) {
	Aux::Timer timer;
	Modularity mod;

	// std::string graph = "../graphs/uk-2002.graph";
	std::string graph = "../graphs/uk-2007-05.graph";
	printf("Reading graph file %s ...\n", graph.c_str());
	timer.start();
	const Graph G = this->metisReader.read(graph);
	timer.stop();
	printf("Reading graph took %.1f s\n", timer.elapsedMilliseconds() / 1000.0);

	std::map<std::string, CommunityDetectionAlgorithm*> algos = {
		std::make_pair("Parallel Label Propagation", (CommunityDetectionAlgorithm*) new PLP),
		std::make_pair("Parallel Louvain", (CommunityDetectionAlgorithm*) new PLM)
	};

	for (auto it = algos.begin(); it != algos.end(); it++) {
		Graph Gcopy = G;

		printf("Timing %s ...\n", it->first.c_str());
		timer.start();
		Partition zeta = it->second->run(Gcopy);
		timer.stop();

		auto communitySizes = zeta.subsetSizes();

		printf("%s on %s: %.1f s\n\t# communities: %lu\n\tmodularity: %f\n",
			it->first.c_str(), graph.c_str(),
			timer.elapsedMilliseconds() / 1000.0,
			zeta.numberOfSubsets(),
			mod.getQuality(zeta, G));
		
		delete it->second;
	}
}

TEST_F(CommunityDetectionBenchmark, timeCentralities) {
	Aux::Timer timer;

	std::string graph = "../graphs/cond-mat-2005.graph";
	const Graph G = this->metisReader.read(graph);

	Graph G1 = this->metisReader.read("../graphs/cond-mat.graph");
	Graph G2 = this->metisReader.read("../graphs/uk-2002.graph");
	std::map<std::string, Centrality*> cens = {
		std::make_pair("Betweenness Centrality on cond-mat", (Centrality*) new Betweenness(G1)),
		std::make_pair("Page Rank Centrality on uk-2002", (Centrality*) new PageRank(G2, 1e-6))
	};

	for (auto it = cens.begin(); it != cens.end(); it++) {
		timer.start();
		it->second->run();
		timer.stop();
		auto ranking = it->second->ranking();

		printf("%s: %.1f s\n\tranking: [(%lu: %f), (%lu: %f), ...]\n",
			it->first.c_str(),
			timer.elapsedMilliseconds() / 1000.0,
			ranking[0].first, ranking[0].second,
			ranking[1].first, ranking[1].second); 
	}
}	

} /* namespace NetworKit */

#endif /*NOGTEST */
