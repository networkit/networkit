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
#include "../../centrality/ApproxBetweenness.h"
#include "../../auxiliary/Timer.h"
#include "../../structures/Partition.h"

namespace NetworKit {

void CommunityDetectionBenchmark::SetUp() {

}

TEST_F(CommunityDetectionBenchmark, timePLP) {
	Aux::Timer timer;
	Modularity mod;
	PLP plp;

	// Graph G = this->metisReader.read("/mnt/windows/Users/Marvin/Desktop/uk-2002.graph");
	Graph G = this->metisReader.read("../in-2004.graph");
	
	timer.start();
	Partition zeta = plp.run(G);
	timer.stop();

	// communitySizes = zeta.subsetSizes();
	printf("PLP: %.2f s\n\t# communities: %lu\n\tmodularity: %f\n",
		timer.elapsedMilliseconds() / 1000.0,
		zeta.numberOfSubsets(),
		mod.getQuality(zeta, G));
}

TEST_F(CommunityDetectionBenchmark, timePLM) {
	Aux::Timer timer;
	Modularity mod;
	PLM plm;

	// Graph G = this->metisReader.read("/mnt/windows/Users/Marvin/Desktop/uk-2002.graph");
	Graph G = this->metisReader.read("../in-2004.graph");

	timer.start();
	Partition zeta = plm.run(G);
	timer.stop();

	// communitySizes = zeta.subsetSizes();
	printf("PLM: %.2f s\n\t# communities: %lu\n\tmodularity: %f\n",
		timer.elapsedMilliseconds() / 1000.0,
		zeta.numberOfSubsets(),
		mod.getQuality(zeta, G));
}

TEST_F(CommunityDetectionBenchmark, timeApproxBetweennessCentrality) {
	Aux::Timer timer;

	Graph G = this->metisReader.read("../cnr-2000.graph");

	timer.start();
	ApproxBetweenness bc(G, 0.1);
	bc.run();
	timer.stop();

	printf("Approximated Betweenness Centrality: %.2f s\n\tscore: [%f, %f, %f, ...]\n",
		timer.elapsedMilliseconds() / 1000.0,
		bc.score(0),
		bc.score(1),
		bc.score(2));
}

} /* namespace NetworKit */

#endif /*NOGTEST */
