/*
 * SpanningGTest.cpp
 *
 *  Created on: 03.09.2015
 *      Author: Henning
 */

#include <gtest/gtest.h>

#include "../KruskalMSF.h"
#include "../SpanningForest.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

class SpanningGTest: public testing::Test {};

TEST_F(SpanningGTest, testKruskalMinSpanningForest) {
	METISGraphReader reader;
	std::vector<std::string> graphs = {"karate", "jazz", "celegans_metabolic"};

	for (auto graphname: graphs) {
		std::string filename = "input/" + graphname + ".graph";
		Graph G = reader.read(filename);
		KruskalMSF msf(G);
		msf.run();
		Graph T = msf.getForest();

		// check that each node has an edge in the spanning tree (if it had one before)
		T.forNodes([&](node u) {
			EXPECT_TRUE(T.degree(u) > 0 || G.degree(u) == 0);
		});
	}
}

TEST_F(SpanningGTest, testSpanningForest) {
	METISGraphReader reader;
	std::vector<std::string> graphs = {"karate", "jazz", "celegans_metabolic"};

	for (auto graphname: graphs) {
		std::string filename = "input/" + graphname + ".graph";
		Graph G = reader.read(filename);
		SpanningForest msf(G);
		Graph T = msf.generate();

		INFO("tree / graph edges: ", T.numberOfEdges(), " / ", G.numberOfEdges());

		// check that each node has an edge in the spanning tree (if it had one before)
		T.forNodes([&](node u) {
//			INFO("tree/graph node degree: ", T.degree(u), " / ", G.degree(u));
			EXPECT_TRUE(T.degree(u) > 0 || G.degree(u) == 0);
		});
	}
}

} /* namespace NetworKit */
