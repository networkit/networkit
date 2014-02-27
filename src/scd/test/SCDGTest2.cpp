#include "SCDGTest2.h"

#include "../PageRankNibble.h"
#include "../../clustering/Modularity.h"
#include "../../clustering/Conductance.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"

#ifndef NOGTEST

namespace NetworKit {

TEST_F(SCDGTest2, testPageRankNibble) {
	METISGraphReader reader;
	Graph G = reader.read("input/polblogs.graph");
	PageRankNibble prn(G);
	count n = G.numberOfNodes();
	count m = G.numberOfEdges();

	// parameters
	node seed = 0;
	double targetCond = 0.4;
	count B = (count) log2(m);
	count b = (count) 0.9 * B;

	// run PageRank-Nibble and partition the graph accordingly
	std::set<node> cluster = prn.run(seed, targetCond, b);
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
