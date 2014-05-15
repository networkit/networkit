#include "SelectiveCDGTest.h"

#include "../PageRankNibble.h"
#include "../SelSCAN.h"
#include "../GCE.h"
#include "../../community/Modularity.h"
#include "../../community/Conductance.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"


#ifndef NOGTEST

namespace NetworKit {

TEST_F(SCDGTest2, testPageRankNibble) {
	METISGraphReader reader;
	Graph G = reader.read("input/hep-th.graph");

	double alpha = 0.1; // loop (or teleport) probability, changed due to DGleich from: // phi * phi / (225.0 * log(100.0 * sqrt(m)));
	double epsilon = 1e-5; // changed due to DGleich from: pow(2, exponent) / (48.0 * B);
	PageRankNibble prn(G, alpha, epsilon);
	count idBound = G.upperNodeIdBound();

	// parameters
	node seed = 50;


	// run PageRank-Nibble and partition the graph accordingly
	DEBUG("Call PageRank-Nibble(", seed, ")");
	std::set<node> cluster = prn.expandSeed({seed});

	// prepare result
	EXPECT_GT(cluster.size(), 0);
	Partition partition(idBound);
	partition.allToOnePartition();
	partition.toSingleton(seed);
	index id = partition[seed];
	for (auto entry: cluster) {
		partition.moveToSubset(id, entry);
	}

	// evaluate result
	Conductance conductance;
	double targetCond = 0.4;
	double cond = conductance.getQuality(partition, G);
	EXPECT_LT(cond, targetCond);
	INFO("Conductance of PR-Nibble: ", cond, "; cluster size: ", cluster.size());
}




TEST_F(SCDGTest2, testSelSCAN) {
	METISGraphReader reader;
	Graph G = reader.read("input/karate.graph");
	count kappa = 2;
	double epsilon = 0.8;

	SelSCAN selSCAN(G, kappa, epsilon);

	auto nodes = G.nodes();
	std::set<unsigned int> seeds(nodes.begin(), nodes.end());

	auto result = selSCAN.run(seeds);

}


TEST_F(SCDGTest2, testGCE) {
	METISGraphReader reader;
	Graph G = reader.read("input/karate.graph");

	auto gce = GCE(G);
	auto nodes = G.nodes();
	std::set<unsigned int> seeds(nodes.begin(), nodes.end());

	auto result = gce.run(seeds);
	DEBUG("communities: ", result);
}


} /* namespace NetworKit */

#endif /*NOGTEST */
