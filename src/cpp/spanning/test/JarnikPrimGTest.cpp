/*
 * JarnikPrimGTest.cpp
 *
 *  Created on: 14.05.2014
 *      Author: Michael
 */

#include "JarnikPrimGTest.h"

namespace NetworKit {

TEST(JarnikPrimGTest, tryMinimumSpanningTree) {
	METISGraphReader graphReader;
	Graph graph = graphReader.read("input/PGPgiantcompo.graph");

	JarnikPrim jp;
	std::vector<std::vector<Edge>> msf = jp.run(graph);
	ASSERT_EQ(1u, msf.size()); // PGPgiantcompo is connected
	ASSERT_EQ(graph.numberOfNodes() - 1, msf.at(0).size());

	std::vector<bool> reachedNodes(graph.numberOfNodes(), false);
	for (auto mst : msf) {
		for (auto edge : mst) {
			reachedNodes.at(edge.first) = true;
			reachedNodes.at(edge.second) = true;
		}
	}

	for (auto reached : reachedNodes) {
		EXPECT_TRUE(reached);
	}
}

TEST(JarnikPrimGTest, tryMinimumSpanningForest) {
	METISGraphReader graphReader;
	Graph graph = graphReader.read("input/polblogs.graph");

	JarnikPrim jp;
	std::vector<std::vector<Edge>> msf = jp.run(graph);
	TRACE("connected components ", msf.size());

	count overallSize = 0;
	for (auto mst : msf) {
		overallSize += mst.size();
	}

	ASSERT_EQ(graph.numberOfNodes() - msf.size(), overallSize);

	std::vector<std::vector<node>> mstNodes;
	std::vector<bool> reachedNodes(graph.numberOfNodes(), false);
	count largestComponent = 0;
	for (auto mst : msf) {
		count componentSize = 0;
		std::vector<node> nodes;
		for (auto edge : mst) {
			if (!reachedNodes.at(edge.first)) {
				componentSize++;
				nodes.push_back(edge.first);
			}

			if (!reachedNodes.at(edge.second)) {
				nodes.push_back(edge.second);
				componentSize++;
			}

			reachedNodes.at(edge.first) = true;
			reachedNodes.at(edge.second) = true;
		}

		mstNodes.push_back(nodes);

		if (componentSize > largestComponent) {
			largestComponent = componentSize;
		}
	}

	TRACE("largest component ", largestComponent);

	// check components
	for (index i = 0; i < mstNodes.size(); ++i) {
		for (auto u : mstNodes.at(i)) {
			graph.forNeighborsOf(u, [&](node v) {
				if (std::find(mstNodes.at(i).begin(), mstNodes.at(i).end(), v) == mstNodes.at(i).end()) {
					TRACE("connected component is incorrect");
				}
			});
		}
	}

	count isolatedNodes = 0;
	graph.forNodes([&](node u){
		if (graph.degree(u) > 0) { // components containing only a single node do not have any edges in mst
			EXPECT_TRUE(reachedNodes[u]);
		} else {
			isolatedNodes++;
		}
	});

	TRACE("isolated nodes ", isolatedNodes);
}


} /* namespace NetworKit */
