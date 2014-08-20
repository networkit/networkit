#include "SelectiveCDGTest.h"

#include "../PageRankNibble.h"
#include "../../community/Modularity.h"
#include "../../community/Conductance.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Log.h"

#ifndef NOGTEST

namespace NetworKit {

TEST_F(SCDGTest2, testPageRankNibble) {
	METISGraphReader reader;
	Graph G = reader.read("input/hep-th.graph");
	PageRankNibble prn(G);
	count idBound = G.upperNodeIdBound();

	// parameters
	node seed = 50;
	double alpha = 0.1; // loop (or teleport) probability, changed due to DGleich from: // phi * phi / (225.0 * log(100.0 * sqrt(m)));
	double epsilon = 1e-5; // changed due to DGleich from: pow(2, exponent) / (48.0 * B);

	// run PageRank-Nibble and partition the graph accordingly
	DEBUG("Call PageRank-Nibble(", seed, ")");
	std::set<node> cluster = prn.run(seed, alpha, epsilon);

	// prepare result
	EXPECT_GT(cluster.size(), 0u);
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


} /* namespace NetworKit */

#endif /*NOGTEST */
