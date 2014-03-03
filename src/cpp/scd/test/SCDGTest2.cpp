#include "SCDGTest2.h"

#include "../PageRankNibble.h"
#include "../../community/Modularity.h"
#include "../../community/Conductance.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"

#ifndef NOGTEST

namespace NetworKit {

TEST_F(SCDGTest2, testPageRankNibble) {
	METISGraphReader reader;
	Graph G = reader.read("input/hep-th.graph");
	PageRankNibble prn(G);
	count n = G.numberOfNodes();
	count m = G.numberOfEdges();

	// parameters
	node seed = 50;
	double targetCond = 0.5;
	count B = (count) ceil(log2(m));
	count b = (count) (0.4 * B); // TODO: vary values!

	// run PageRank-Nibble and partition the graph accordingly
	DEBUG("Call PageRank-Nibble(", seed, ", ", targetCond, ", ", b, "), B=", B);
	std::set<node> cluster = prn.run(seed, targetCond, b);
	EXPECT_GT(cluster.size(), 0);
	Partition partition(n, 0);
	for (auto entry: cluster) {
		partition[entry] = 1;
	}

	// evaluate result
	Modularity modularity;
	Conductance conductance;
	double mod = modularity.getQuality(partition, G);
	double cond = conductance.getQuality(partition, G);

	EXPECT_LT(cond, targetCond);
	INFO("Conductance: ", cond);
	INFO("Modularity: ", mod);
}


} /* namespace NetworKit */

#endif /*NOGTEST */
