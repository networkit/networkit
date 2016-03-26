#include "MaximalCliquesGTest.h"

#include "../MaximalCliques.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Log.h"

#ifndef NOGTEST

namespace NetworKit {

TEST_F(MaximalCliquesGTest, testMaximalCliques) {

	METISGraphReader reader;
	Graph G = reader.read("input/hep-th.graph");

	node seed = 2;
	auto sn = G.neighbors(seed);
	auto sneighbors = std::unordered_set<node>(sn.begin(), sn.end());
	sneighbors.insert(seed);
	auto subG = G.subgraphFromNodes(sneighbors);

	MaximalCliques clique(subG);

	DEBUG("Call MaximalCliques()");

	auto result = clique.run();

	EXPECT_GT(result.size(), 0);

	// check results (are they cliques?)
	for (auto cliq : result) {
		auto cli = std::unordered_set<node>(cliq.begin(), cliq.end());
		auto cliqueGraph = Subgraph::fromNodes(G, cli);

		EXPECT_EQ(cliqueGraph.numberOfEdges(), (cliqueGraph.numberOfNodes() * (cliqueGraph.numberOfNodes() - 1) / 2));
		EXPECT_EQ(cli.count(seed), 1);
	}
}

}

#endif /* NOGTEST */
