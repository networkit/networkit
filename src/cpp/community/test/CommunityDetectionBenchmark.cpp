/*
 * CommunityDetectionBenchmark.h
 *
 *  Created on: 16.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

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

TEST_F(CommunityDetectionBenchmark, timePLP) {
	Aux::Timer timer;

	std::string graph = "~/graphs/in-2004.graph";
	Graph G = this->metisReader.read(graph);

	std::string name = "PLP";
	CommunityDetectionAlgorithm* algo = new PLP;

	timer.start();
	Partition zeta = algo->run(G);
	timer.stop();

	Modularity mod;
	auto communitySizes = zeta.subsetSizes();
	delete algo;

	printf("%s on %s: %.1f s\n\t# communities: %lu\n\tmodularity: %f\n",
		name.c_str(), graph.c_str(),
		timer.elapsedMilliseconds() / 1000.0,
		zeta.numberOfSubsets(),
		mod.getQuality(zeta, G));
}

TEST_F(CommunityDetectionBenchmark, timePLM) {

	std::string graph = "~/graphs/in-2004.graph";
	Graph G = this->metisReader.read(graph);

	std::string name = "PLM";
	CommunityDetectionAlgorithm* algo = new PLM;

	timer.start();
	Partition zeta = algo->run(G);
	timer.stop();

	Modularity mod;
	auto communitySizes = zeta.subsetSizes();
	delete algo;

	printf("%s on %s: %.1f s\n\t# communities: %lu\n\tmodularity: %f\n",
		name.c_str(), graph.c_str(),
		timer.elapsedMilliseconds() / 1000.0,
		zeta.numberOfSubsets(),
		mod.getQuality(zeta, G));
}

TEST_F(CommunityDetectionBenchmark, timeBetweennessCentrality) {
	Aux::Timer timer;


	std::string graph = "~/graphs/in-2004.graph";
	Graph G = this->metisReader.read(graph);

	Centrality* cen = new Betweenness(G, 1e-6);
	std::string name = "Betweenness Centrality";

	timer.start();
	cen->run();
	timer.stop();
	auto ranking = cen->ranking();
	delete cen;

	printf("%s on %s: %.1f s\n\tranking: [(%lu: %f), (%lu: %f), ...]\n",
		name.c_str(), graph.c_str(),
		timer.elapsedMilliseconds() / 1000.0,
		ranking[0].first, ranking[0].second,
		ranking[1].first, ranking[1].second); 
}

TEST_F(CommunityDetectionBenchmark, timePageRankCentrality) {
	Aux::Timer timer;

	std::string graph = "~/graphs/in-2004.graph";
	Graph G = this->metisReader.read(graph);

	Centrality* cen = new PageRank(G, 1e-6);
	std::string name = "Page Rank Centrality";

	timer.start();
	cen->run();
	timer.stop();
	auto ranking = cen->ranking();
	delete cen;

	printf("%s on %s: %.1f s\n\tranking: [(%lu: %f), (%lu: %f), ...]\n",
		name.c_str(), graph.c_str(),
		timer.elapsedMilliseconds() / 1000.0,
		ranking[0].first, ranking[0].second,
		ranking[1].first, ranking[1].second); 
}

} /* namespace NetworKit */

#endif /*NOGTEST */
