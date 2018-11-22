/*
 * GraphGTest.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter
 * (marvin.ritter@gmail.com)
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <tuple>

#include "../Graph.h"

#include "../../auxiliary/NumericTools.h"
#include "../../auxiliary/Parallel.h"
#include "../../distance/DynBFS.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../../io/METISGraphReader.h"
#include "../GraphBuilder.h"

namespace NetworKit {

class GraphGTest : public testing::TestWithParam<std::tuple<bool, bool>> {
public:
	virtual void SetUp();

protected:
	Graph Ghouse;
	std::vector<std::pair<node, node>> houseEdgesOut;
	std::vector<std::vector<edgeweight>> Ahouse;
	count n_house;
	count m_house;

	bool isGraph() const { return !isWeighted() && !isDirected(); }
	bool isWeightedGraph() const { return isWeighted() && !isDirected(); }
	bool isDirectedGraph() const { return !isWeighted() && isDirected(); }
	bool isWeightedDirectedGraph() const { return isWeighted() && isDirected(); }

	bool isWeighted() const;
	bool isDirected() const;
	Graph createGraph(count n = 0) const;
	count countSelfLoopsManually(const Graph &G);
};

INSTANTIATE_TEST_CASE_P(InstantiationName, GraphGTest,
                        testing::Values(std::make_tuple(false, false),
                                        std::make_tuple(true, false),
                                        std::make_tuple(false, true),
                                        std::make_tuple(true, true)));

bool GraphGTest::isWeighted() const { return std::get<0>(GetParam()); }
bool GraphGTest::isDirected() const { return std::get<1>(GetParam()); }

Graph GraphGTest::createGraph(count n) const {
	bool weighted, directed;
	std::tie(weighted, directed) = GetParam();
	Graph G(n, weighted, directed);
	return G;
}

count GraphGTest::countSelfLoopsManually(const Graph &G) {
	count c = 0;
	G.parallelForEdges([&](node u, node v) {
		if (u == v) {
#pragma omp atomic
			c += 1;
		}
	});
	return c;
}

void GraphGTest::SetUp() {
	/*
	 *    0
	 *   . \
	 *  /   \
	 * /     .
	 * 1 <-- 2
	 * ^ \  .|
	 * |  \/ |
	 * | / \ |
	 * |/   ..
	 * 3 <-- 4
	 *
	 * move you pen from node to node:
	 * 3 -> 1 -> 0 -> 2 -> 1 -> 4 -> 3 -> 2 -> 4
	 */
	n_house = 5;
	m_house = 8;

	Ghouse = createGraph(5);
	houseEdgesOut = {{3, 1}, {1, 0}, {0, 2}, {2, 1},
	                 {1, 4}, {4, 3}, {3, 2}, {2, 4}};
	Ahouse = {n_house, std::vector<edgeweight>(n_house, 0.0)};
	edgeweight ew = 1.0;
	for (auto &e : houseEdgesOut) {
		node u = e.first;
		node v = e.second;
		Ghouse.addEdge(u, v, ew);

		Ahouse[u][v] = ew;

		if (!Ghouse.isDirected()) {
			Ahouse[v][u] = ew;
		}

		if (Ghouse.isWeighted()) {
			ew += 1.0;
		}
	}
}

/** CONSTRUCTORS **/

TEST_P(GraphGTest, testCopyConstructor) {
	Graph G = Graph(this->Ghouse, false, false);
	Graph GW = Graph(this->Ghouse, true, false);
	Graph D = Graph(this->Ghouse, false, true);
	Graph DW = Graph(this->Ghouse, true, true);

	ASSERT_FALSE(G.isWeighted());
	ASSERT_FALSE(G.isDirected());
	ASSERT_EQ(this->Ghouse.numberOfNodes(), G.numberOfNodes());
	ASSERT_EQ(this->Ghouse.numberOfEdges(), G.numberOfEdges());

	ASSERT_TRUE(GW.isWeighted());
	ASSERT_FALSE(GW.isDirected());
	ASSERT_EQ(this->Ghouse.numberOfNodes(), GW.numberOfNodes());
	ASSERT_EQ(this->Ghouse.numberOfEdges(), GW.numberOfEdges());

	ASSERT_FALSE(D.isWeighted());
	ASSERT_TRUE(D.isDirected());
	ASSERT_EQ(this->Ghouse.numberOfNodes(), D.numberOfNodes());
	ASSERT_EQ(this->Ghouse.numberOfEdges(), D.numberOfEdges());

	ASSERT_TRUE(DW.isWeighted());
	ASSERT_TRUE(DW.isDirected());
	ASSERT_EQ(this->Ghouse.numberOfNodes(), DW.numberOfNodes());
	ASSERT_EQ(this->Ghouse.numberOfEdges(), DW.numberOfEdges());

	this->Ghouse.forNodes([&](node v) {
		count d = this->Ghouse.degree(v);
		count dUndirected = isDirected() ? d + this->Ghouse.degreeIn(v) : d;
		ASSERT_EQ(dUndirected, G.degree(v));
		ASSERT_EQ(dUndirected, GW.degree(v));
		ASSERT_EQ(d, D.degree(v));
		ASSERT_EQ(d, DW.degree(v));
	});

	// if Ghouse was directed we should have an exact copy of it, but if it was
	// undirected we should have edges in both directions
	count m = 0;
	G.forEdges([&](node u, node v) {
		ASSERT_TRUE(G.hasEdge(v, u));
		ASSERT_EQ(defaultEdgeWeight, G.weight(v, u));
		ASSERT_EQ(defaultEdgeWeight, G.weight(u, v));

		auto e = std::make_pair(u, v);
		bool found =
		    std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) !=
		    this->houseEdgesOut.end();
		if (!found) {
			e = std::make_pair(v, u);
			found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
			                  e) != this->houseEdgesOut.end();
		}
		ASSERT_TRUE(found);
		m++;
	});
	ASSERT_EQ(8u, m);

	m = 0;
	GW.forEdges([&](node u, node v) {
		ASSERT_TRUE(GW.hasEdge(v, u));
		ASSERT_EQ(GW.weight(u, v), GW.weight(v, u));

		auto e = std::make_pair(u, v);
		bool found =
		    std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) !=
		    this->houseEdgesOut.end();
		if (!found) {
			e = std::make_pair(v, u);
			found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
			                  e) != this->houseEdgesOut.end();
			ASSERT_EQ(this->Ghouse.weight(v, u), GW.weight(v, u));
		} else {
			ASSERT_EQ(this->Ghouse.weight(u, v), GW.weight(u, v));
		}
		ASSERT_TRUE(found);
		m++;
	});
	ASSERT_EQ(8u, m);

	m = 0;
	D.forEdges([&](node u, node v) {
		ASSERT_EQ(defaultEdgeWeight, D.weight(u, v));

		auto e = std::make_pair(u, v);
		bool found =
		    std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) !=
		    this->houseEdgesOut.end();
		if (!this->Ghouse.isDirected()) {
			ASSERT_TRUE(D.hasEdge(v, u));

			e = std::make_pair(v, u);
			found = found ||
			        (std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
			                   e) != this->houseEdgesOut.end());
		} else {
			ASSERT_FALSE(D.hasEdge(v, u));
		}
		ASSERT_TRUE(found);
		m++;
	});
	count m_expected = isDirected() ? 8 : 16;
	ASSERT_EQ(m_expected, m);

	m = 0;
	DW.forEdges([&](node u, node v) {
		ASSERT_EQ(this->Ghouse.weight(u, v), DW.weight(u, v));

		auto e = std::make_pair(u, v);
		bool found =
		    std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) !=
		    this->houseEdgesOut.end();
		if (!this->Ghouse.isDirected()) {
			ASSERT_TRUE(DW.hasEdge(v, u));
			e = std::make_pair(v, u);
			found = found ||
			        (std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
			                   e) != this->houseEdgesOut.end());
		} else {
			ASSERT_FALSE(DW.hasEdge(v, u));
		}
		ASSERT_TRUE(found);
		m++;
	});
	m_expected = isDirected() ? 8 : 16;
	ASSERT_EQ(m_expected, m);
}

/** GRAPH INFORMATION **/

TEST_P(GraphGTest, testGetId) {
	Graph G1 = createGraph();
	Graph G2 = createGraph(5);

	ASSERT_TRUE(G1.getId() > 0);
	ASSERT_TRUE(G2.getId() > 0);
	ASSERT_TRUE(G1.getId() < G2.getId());
}

TEST_P(GraphGTest, testTyp) {
	Graph G = createGraph();
	if (isGraph()) {
		ASSERT_EQ("Graph", G.typ());
	} else if (isWeightedGraph()) {
		ASSERT_EQ("WeightedGraph", G.typ());
	} else if (isDirectedGraph()) {
		ASSERT_EQ("DirectedGraph", G.typ());
	} else if (isWeightedDirectedGraph()) {
		ASSERT_EQ("WeightedDirectedGraph", G.typ());
	} else {
		FAIL();
	}
}

TEST_P(GraphGTest, testSetName) {
	Graph G1 = createGraph(0);
	Graph G2 = createGraph(0);

	std::string s1 = "Graph 1";
	std::string s2 = "Graph 2";
	G1.setName(s1);
	G2.setName(s2);
	ASSERT_EQ(s1, G1.getName());
	ASSERT_EQ(s2, G2.getName());
}

TEST_P(GraphGTest, testMaxDegreeUndirected) {
	Aux::Random::setSeed(1, false);
	Graph G = ErdosRenyiGenerator(20, 0.2, false).generate();

	count maxDegOut = 0, maxDegIn = 0;
	G.forNodes([&](const node u) {
		maxDegOut = std::max(maxDegOut, G.degreeOut(u));
		maxDegIn = std::max(maxDegIn, G.degreeIn(u));
	});

	ASSERT_EQ(G.maxDegree(), maxDegOut);
	ASSERT_EQ(G.maxDegreeIn(), maxDegIn);
}

TEST_P(GraphGTest, testMaxDegreeDirected) {
	Aux::Random::setSeed(1, false);
	Graph G = ErdosRenyiGenerator(20, 0.2, true).generate();

	count maxDegOut = 0, maxDegIn = 0;
	G.forNodes([&](const node u) {
		maxDegOut = std::max(maxDegOut, G.degreeOut(u));
		maxDegIn = std::max(maxDegIn, G.degreeIn(u));
	});

	ASSERT_EQ(G.maxDegree(), maxDegOut);
	ASSERT_EQ(G.maxDegreeIn(), maxDegIn);
}

TEST_P(GraphGTest, testMaxWeightedDegreeUndirected) {
	Aux::Random::setSeed(1, false);
	Graph G = ErdosRenyiGenerator(20, 0.2, false).generate();

	edgeweight maxDegOut = 0.0, maxDegIn = 0.0;
	G.forNodes([&](const node u) {
		maxDegOut = std::max(maxDegOut, G.weightedDegree(u));
		maxDegIn = std::max(maxDegIn, G.weightedDegreeIn(u));
	});

	ASSERT_EQ(G.maxWeightedDegree(), maxDegOut);
	ASSERT_EQ(G.maxWeightedDegreeIn(), maxDegIn);
}

TEST_P(GraphGTest, testMaxWeightedDegreeDirected) {
	Aux::Random::setSeed(1, false);
	Graph G = ErdosRenyiGenerator(20, 0.2, true).generate();

	edgeweight maxDegOut = 0.0, maxDegIn = 0.0;
	G.forNodes([&](const node u) {
		maxDegOut = std::max(maxDegOut, G.weightedDegree(u));
		maxDegIn = std::max(maxDegIn, G.weightedDegreeIn(u));
	});

	ASSERT_EQ(G.maxWeightedDegree(), maxDegOut);
	ASSERT_EQ(G.maxWeightedDegreeIn(), maxDegIn);
}

TEST_P(GraphGTest, testToString) {
	Graph G1 = createGraph(0);
	Graph G2 = createGraph(0);

	ASSERT_TRUE(G1.toString() != "");
	ASSERT_TRUE(G2.toString() != "");
}

/** NODE MODIFIERS **/

TEST_P(GraphGTest, testAddNode) {
	Graph G = createGraph();

	ASSERT_FALSE(G.hasNode(0));
	ASSERT_FALSE(G.hasNode(1));
	ASSERT_EQ(0u, G.numberOfNodes());

	G.addNode();
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_FALSE(G.hasNode(1));
	ASSERT_EQ(1u, G.numberOfNodes());

	Graph G2 = createGraph(2);
	ASSERT_TRUE(G2.hasNode(0));
	ASSERT_TRUE(G2.hasNode(1));
	ASSERT_FALSE(G2.hasNode(2));
	ASSERT_EQ(2u, G2.numberOfNodes());

	G2.addNode();
	G2.addNode();
	ASSERT_TRUE(G2.hasNode(2));
	ASSERT_TRUE(G2.hasNode(3));
	ASSERT_FALSE(G2.hasNode(4));
	ASSERT_EQ(4u, G2.numberOfNodes());
}

TEST_P(GraphGTest, testRemoveNode) {
	Graph G = createGraph(4);
	G.addEdge(0, 1);
	G.addEdge(0, 2);

	ASSERT_EQ(4u, G.numberOfNodes());
	ASSERT_EQ(2u, G.numberOfEdges());
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_TRUE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));

	G.removeNode(2);

	ASSERT_EQ(3u, G.numberOfNodes());
	ASSERT_EQ(1u, G.numberOfEdges());
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_FALSE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));

	G.removeNode(1);

	ASSERT_EQ(2u, G.numberOfNodes());
	ASSERT_EQ(0u, G.numberOfEdges());
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_FALSE(G.hasNode(1));
	ASSERT_FALSE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));
}

TEST_P(GraphGTest, testHasNode) {
	Graph G = createGraph(5);

	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_TRUE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));
	ASSERT_TRUE(G.hasNode(4));
	ASSERT_FALSE(G.hasNode(5));
	ASSERT_FALSE(G.hasNode(6));

	G.removeNode(0);
	G.removeNode(2);
	G.addNode();

	ASSERT_FALSE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_FALSE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));
	ASSERT_TRUE(G.hasNode(4));
	ASSERT_TRUE(G.hasNode(5));
	ASSERT_FALSE(G.hasNode(6));
}

TEST_P(GraphGTest, testRestoreNode) {
	Graph G = createGraph(4);

	ASSERT_EQ(4u, G.numberOfNodes());
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_TRUE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));

	G.removeNode(0);

	ASSERT_EQ(3u, G.numberOfNodes());
	ASSERT_FALSE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_TRUE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));

	G.restoreNode(0);

	ASSERT_EQ(4u, G.numberOfNodes());
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_TRUE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));
}

/** NODE PROPERTIES **/

TEST_P(GraphGTest, testDegree) {
	if (isDirected()) {
		ASSERT_EQ(1u, this->Ghouse.degree(0));
		ASSERT_EQ(2u, this->Ghouse.degree(1));
		ASSERT_EQ(2u, this->Ghouse.degree(2));
		ASSERT_EQ(2u, this->Ghouse.degree(3));
		ASSERT_EQ(1u, this->Ghouse.degree(4));
	} else {
		ASSERT_EQ(2u, this->Ghouse.degree(0));
		ASSERT_EQ(4u, this->Ghouse.degree(1));
		ASSERT_EQ(4u, this->Ghouse.degree(2));
		ASSERT_EQ(3u, this->Ghouse.degree(3));
		ASSERT_EQ(3u, this->Ghouse.degree(4));
	}
}

TEST_P(GraphGTest, testDegreeIn) {
	if (isDirected()) {
		ASSERT_EQ(1u, this->Ghouse.degreeIn(0));
		ASSERT_EQ(2u, this->Ghouse.degreeIn(1));
		ASSERT_EQ(2u, this->Ghouse.degreeIn(2));
		ASSERT_EQ(1u, this->Ghouse.degreeIn(3));
		ASSERT_EQ(2u, this->Ghouse.degreeIn(4));
	} else {
		ASSERT_EQ(2u, this->Ghouse.degreeIn(0));
		ASSERT_EQ(4u, this->Ghouse.degreeIn(1));
		ASSERT_EQ(4u, this->Ghouse.degreeIn(2));
		ASSERT_EQ(3u, this->Ghouse.degreeIn(3));
		ASSERT_EQ(3u, this->Ghouse.degreeIn(4));
	}
}

TEST_P(GraphGTest, testDegreeOut) {
	if (isDirected()) {
		ASSERT_EQ(1u, this->Ghouse.degreeOut(0));
		ASSERT_EQ(2u, this->Ghouse.degreeOut(1));
		ASSERT_EQ(2u, this->Ghouse.degreeOut(2));
		ASSERT_EQ(2u, this->Ghouse.degreeOut(3));
		ASSERT_EQ(1u, this->Ghouse.degreeOut(4));
	} else {
		ASSERT_EQ(2u, this->Ghouse.degreeOut(0));
		ASSERT_EQ(4u, this->Ghouse.degreeOut(1));
		ASSERT_EQ(4u, this->Ghouse.degreeOut(2));
		ASSERT_EQ(3u, this->Ghouse.degreeOut(3));
		ASSERT_EQ(3u, this->Ghouse.degreeOut(4));
	}
}

TEST_P(GraphGTest, testIsIsolated) {
	ASSERT_FALSE(this->Ghouse.isIsolated(0));
	ASSERT_FALSE(this->Ghouse.isIsolated(1));
	ASSERT_FALSE(this->Ghouse.isIsolated(2));
	ASSERT_FALSE(this->Ghouse.isIsolated(3));
	ASSERT_FALSE(this->Ghouse.isIsolated(4));

	this->Ghouse.addNode();
	ASSERT_TRUE(this->Ghouse.isIsolated(5));

	this->Ghouse.removeEdge(1, 0);
	ASSERT_FALSE(this->Ghouse.isIsolated(0));

	this->Ghouse.removeEdge(0, 2);
	ASSERT_TRUE(this->Ghouse.isIsolated(0));

	this->Ghouse.addEdge(1, 0);
	ASSERT_FALSE(this->Ghouse.isIsolated(0));
}

TEST_P(GraphGTest, testWeightedDegree) {
	// add self-loop
	this->Ghouse.addEdge(2, 2, 0.75);

	if (isGraph()) {
		ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(0));
		ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.weightedDegree(1));
		ASSERT_EQ(5 * defaultEdgeWeight, this->Ghouse.weightedDegree(2));
		ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(3));
		ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(4));
	}

	if (isWeightedGraph()) {
		ASSERT_EQ(5.0, this->Ghouse.weightedDegree(0));
		ASSERT_EQ(12.0, this->Ghouse.weightedDegree(1));
		ASSERT_EQ(22.75, this->Ghouse.weightedDegree(2));
		ASSERT_EQ(14.0, this->Ghouse.weightedDegree(3));
		ASSERT_EQ(19.0, this->Ghouse.weightedDegree(4));
	}

	if (isDirectedGraph()) {
		// only count outgoing edges
		ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(0));
		ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(1));
		ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(2));
		ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(3));
		ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(4));
	}

	if (isWeightedDirectedGraph()) {
		// only sum weight of outgoing edges
		ASSERT_EQ(3.0, this->Ghouse.weightedDegree(0));
		ASSERT_EQ(7.0, this->Ghouse.weightedDegree(1));
		ASSERT_EQ(12.75, this->Ghouse.weightedDegree(2));
		ASSERT_EQ(8.0, this->Ghouse.weightedDegree(3));
		ASSERT_EQ(6.0, this->Ghouse.weightedDegree(4));
	}
}

TEST_P(GraphGTest, testVolume) {
	// add self-loop
	this->Ghouse.addEdge(2, 2, 0.75);

	if (isGraph()) {
		ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.volume(0));
		ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.volume(1));
		ASSERT_EQ(6 * defaultEdgeWeight, this->Ghouse.volume(2));
		ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.volume(3));
		ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.volume(4));
	}

	if (isWeightedGraph()) {
		ASSERT_EQ(5.0, this->Ghouse.volume(0));
		ASSERT_EQ(12.0, this->Ghouse.volume(1));
		ASSERT_EQ(23.5, this->Ghouse.volume(2));
		ASSERT_EQ(14.0, this->Ghouse.volume(3));
		ASSERT_EQ(19.0, this->Ghouse.volume(4));
	}

	if (isDirectedGraph()) {
		// only count outgoing edges
		ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.volume(0));
		ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.volume(1));
		ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.volume(2));
		ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.volume(3));
		ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.volume(4));
	}

	if (isWeightedDirectedGraph()) {
		// only sum weight of outgoing edges
		ASSERT_EQ(3.0, this->Ghouse.volume(0));
		ASSERT_EQ(7.0, this->Ghouse.volume(1));
		ASSERT_EQ(13.5, this->Ghouse.volume(2));
		ASSERT_EQ(8.0, this->Ghouse.volume(3));
		ASSERT_EQ(6.0, this->Ghouse.volume(4));
	}
}

TEST_P(GraphGTest, testRandomNode) {
	count n = 4;
	count samples = 100000;
	double maxAbsoluteError = 0.005;
	Aux::Random::setSeed(42, false);

	Graph G = createGraph(n);
	std::vector<count> drawCounts(n, 0);
	for (count i = 0; i < samples; i++) {
		node x = G.randomNode();
		drawCounts[x]++;
	}
	for (node v = 0; v < n; v++) {
		double p = drawCounts[v] / (double)samples;
		ASSERT_NEAR(1.0 / n, p, maxAbsoluteError);
	}
}

TEST_P(GraphGTest, testRandomNeighbor) {
	Graph G = createGraph(10);
	G.addEdge(2, 0);
	G.addEdge(2, 1);
	G.addEdge(2, 2);
	G.addEdge(5, 6);

	Aux::Random::setSeed(42, false);

	ASSERT_EQ(none, G.randomNeighbor(3));
	ASSERT_EQ(6u, G.randomNeighbor(5));

	if (G.isDirected()) {
		ASSERT_EQ(none, G.randomNeighbor(1));
	} else {
		ASSERT_EQ(2u, G.randomNeighbor(1));
	}

	count nn = 3;
	count samples = 100000;
	double maxAbsoluteError = 0.005;
	std::vector<count> drawCounts(nn, 0);
	for (count i = 0; i < samples; i++) {
		node x = G.randomNeighbor(2);
		drawCounts[x]++;
	}
	for (node v = 0; v < nn; v++) {
		double p = drawCounts[v] / (double)samples;
		ASSERT_NEAR(1.0 / nn, p, maxAbsoluteError);
	}
}

/** EDGE MODIFIERS **/

TEST_P(GraphGTest, testAddEdge) {
	Graph G = createGraph(3);

	// Graph without edges
	ASSERT_EQ(0u, G.numberOfEdges());
	ASSERT_FALSE(G.hasEdge(0, 2));
	ASSERT_FALSE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(1, 2));
	ASSERT_FALSE(G.hasEdge(2, 2));
	ASSERT_EQ(nullWeight, G.weight(0, 2));
	ASSERT_EQ(nullWeight, G.weight(0, 1));
	ASSERT_EQ(nullWeight, G.weight(1, 2));
	ASSERT_EQ(nullWeight, G.weight(2, 2));

	// Graph with 2 normal edges
	G.addEdge(0, 1, 4.51);
	G.addEdge(1, 2, 2.39);
	ASSERT_EQ(2u, G.numberOfEdges());
	ASSERT_FALSE(G.hasEdge(0, 2)); // was never added
	ASSERT_TRUE(G.hasEdge(0, 1));
	ASSERT_TRUE(G.hasEdge(1, 2));
	ASSERT_FALSE(G.hasEdge(2, 2)); // will be added later

	// check weights
	if (G.isWeighted()) {
		ASSERT_EQ(4.51, G.weight(0, 1));
		ASSERT_EQ(2.39, G.weight(1, 2));
	} else {
		ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
		ASSERT_EQ(defaultEdgeWeight, G.weight(1, 2));
	}

	if (G.isDirected()) {
		ASSERT_FALSE(G.hasEdge(1, 0));
		ASSERT_FALSE(G.hasEdge(2, 1));

		// add edge in the other direction
		// note: bidirectional edges are not supported, so both edges have different
		// weights
		G.addEdge(2, 1, 6.23);
		ASSERT_TRUE(G.hasEdge(2, 1));
		if (G.isWeighted()) {
			ASSERT_EQ(2.39, G.weight(1, 2));
			ASSERT_EQ(6.23, G.weight(2, 1));
		} else {
			ASSERT_EQ(defaultEdgeWeight, G.weight(2, 1));
		}
	} else {
		ASSERT_TRUE(G.hasEdge(1, 0));
		ASSERT_TRUE(G.hasEdge(2, 1));
		if (G.isWeighted()) {
			ASSERT_EQ(4.51, G.weight(1, 0));
			ASSERT_EQ(2.39, G.weight(2, 1));
		} else {
			ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
			ASSERT_EQ(defaultEdgeWeight, G.weight(2, 1));
		}
	}

	// add self loop
	G.addEdge(2, 2, 0.72);
	ASSERT_TRUE(G.hasEdge(2, 2));
	if (G.isWeighted()) {
		ASSERT_EQ(0.72, G.weight(2, 2));
	} else {
		ASSERT_EQ(defaultEdgeWeight, G.weight(2, 2));
	}
}

TEST_P(GraphGTest, testRemoveEdge) {
	double epsilon = 1e-6;
	Graph G = createGraph(3);

	edgeweight ewBefore = G.totalEdgeWeight();

	G.addEdge(0, 1, 3.14);

	if (G.isWeighted()) {
		ASSERT_NEAR(ewBefore + 3.14, G.totalEdgeWeight(), epsilon);
	} else {
		ASSERT_NEAR(ewBefore + defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
	}

	G.addEdge(0, 0);

	ASSERT_EQ(2u, G.numberOfEdges());
	ASSERT_TRUE(G.hasEdge(0, 0));
	ASSERT_TRUE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));

	// test remove regular edge
	ewBefore = G.totalEdgeWeight();
	G.removeEdge(0, 1);
	if (G.isWeighted()) {
		ASSERT_NEAR(ewBefore - 3.14, G.totalEdgeWeight(), epsilon);
	} else {
		ASSERT_NEAR(ewBefore - defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
	}

	ASSERT_EQ(1u, G.numberOfEdges());
	ASSERT_TRUE(G.hasEdge(0, 0));
	ASSERT_FALSE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));

	// test remove self-loop
	G.addEdge(2, 1);

	ewBefore = G.totalEdgeWeight();
	G.removeEdge(0, 0);
	if (G.isWeighted()) {
		ASSERT_NEAR(ewBefore - defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
	} else {
		ASSERT_NEAR(ewBefore - defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
	}

	ASSERT_EQ(1u, G.numberOfEdges());
	ASSERT_FALSE(G.hasEdge(0, 0));
	ASSERT_FALSE(G.hasEdge(0, 1));
	ASSERT_TRUE(G.hasEdge(2, 1));

	// test from removeselfloops adapted for removeEdge
	G = createGraph(2);

	ewBefore = G.totalEdgeWeight();

	G.addEdge(0, 1);
	G.addEdge(0, 0, 3.14);
	G.addEdge(1, 1);

	if (G.isWeighted()) {
		EXPECT_NEAR(ewBefore + 3.14 + 2 * defaultEdgeWeight, G.totalEdgeWeight(),
		            epsilon);
	} else {
		EXPECT_NEAR(ewBefore + 3 * defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
	}

	EXPECT_EQ(3u, G.numberOfEdges());
	EXPECT_TRUE(G.hasEdge(0, 0));
	EXPECT_TRUE(G.hasEdge(0, 1));
	EXPECT_TRUE(G.hasEdge(1, 1));
	EXPECT_EQ(G.numberOfSelfLoops(), 2u);

	// remove self-loops
	ewBefore = G.totalEdgeWeight();

	// G.removeSelfLoops();
	G.removeEdge(0, 0);
	G.removeEdge(1, 1);

	if (G.isWeighted()) {
		EXPECT_NEAR(ewBefore - defaultEdgeWeight - 3.14, G.totalEdgeWeight(),
		            epsilon);
	} else {
		EXPECT_NEAR(ewBefore - 2 * defaultEdgeWeight, G.totalEdgeWeight(), epsilon)
		    << "Weighted, directed: " << G.isWeighted() << ", " << G.isDirected();
	}

	EXPECT_EQ(1u, G.numberOfEdges());
	EXPECT_FALSE(G.hasEdge(0, 0));
	EXPECT_FALSE(G.hasEdge(1, 1));
	EXPECT_TRUE(G.hasEdge(0, 1));
	EXPECT_EQ(0u, G.numberOfSelfLoops())
	    << "Weighted, directed: " << G.isWeighted() << ", " << G.isDirected();
}

TEST_P(GraphGTest, testRemoveAllEdges) {
	Graph g = ErdosRenyiGenerator(20, 0.1, false).generate();
	g.removeAllEdges();
	EXPECT_EQ(g.numberOfEdges(), 0);
	EXPECT_EQ(g.edges().size(), 0);
	for (node u : g.nodes()) {
		EXPECT_EQ(g.neighbors(u).size(), 0);
		EXPECT_EQ(g.degree(u), 0);
	}

	g = ErdosRenyiGenerator(20, 0.1, true).generate();
	g.removeAllEdges();
	EXPECT_EQ(g.numberOfEdges(), 0);
	EXPECT_EQ(g.edges().size(), 0);
	for (node u : g.nodes()) {
		EXPECT_EQ(g.neighbors(u).size(), 0);
		EXPECT_EQ(g.degree(u), 0);
		EXPECT_EQ(g.degreeIn(u), 0);
	}
}

TEST_P(GraphGTest, testRemoveSelfLoops) {
	double epsilon = 1e-6;
	Graph G = createGraph(2);

	edgeweight ewBefore = G.totalEdgeWeight();

	G.addEdge(0, 1);
	G.addEdge(0, 0, 3.14);
	G.addEdge(1, 1);

	if (G.isWeighted()) {
		EXPECT_NEAR(ewBefore + 3.14 + 2 * defaultEdgeWeight, G.totalEdgeWeight(),
		            epsilon);
	} else {
		EXPECT_NEAR(ewBefore + 3 * defaultEdgeWeight, G.totalEdgeWeight(), epsilon);
	}

	EXPECT_EQ(3u, G.numberOfEdges());
	EXPECT_TRUE(G.hasEdge(0, 0));
	EXPECT_TRUE(G.hasEdge(0, 1));
	EXPECT_TRUE(G.hasEdge(1, 1));
	EXPECT_EQ(G.numberOfSelfLoops(), 2u);
	EXPECT_EQ(G.numberOfSelfLoops(), 2u);

	// remove self-loops
	ewBefore = G.totalEdgeWeight();
	G.removeSelfLoops();

	if (G.isWeighted()) {
		EXPECT_NEAR(ewBefore - defaultEdgeWeight - 3.14, G.totalEdgeWeight(),
		            epsilon);
	} else {
		EXPECT_NEAR(ewBefore - 2 * defaultEdgeWeight, G.totalEdgeWeight(), epsilon)
		    << "Weighted, directed: " << G.isWeighted() << ", " << G.isDirected();
	}

	EXPECT_EQ(1u, G.numberOfEdges());
	EXPECT_FALSE(G.hasEdge(0, 0));
	EXPECT_FALSE(G.hasEdge(1, 1));
	EXPECT_TRUE(G.hasEdge(0, 1));
	EXPECT_EQ(0u, G.numberOfSelfLoops())
	    << "Weighted, directed: " << G.isWeighted() << ", " << G.isDirected();
}

TEST_P(GraphGTest, testHasEdge) {
	auto containsEdge = [&](std::pair<node, node> e) {
		auto it =
		    std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
		return it != this->houseEdgesOut.end();
	};

	for (node u = 0; u < this->Ghouse.upperNodeIdBound(); u++) {
		for (node v = 0; v < this->Ghouse.upperNodeIdBound(); v++) {
			auto edge = std::make_pair(u, v);
			auto edgeReverse = std::make_pair(v, u);
			bool hasEdge = containsEdge(edge);
			bool hasEdgeReverse = containsEdge(edgeReverse);
			if (this->Ghouse.isDirected()) {
				ASSERT_EQ(hasEdge, this->Ghouse.hasEdge(u, v));
			} else {
				ASSERT_EQ(hasEdge || hasEdgeReverse, this->Ghouse.hasEdge(u, v));
			}
		}
	}
}

TEST_P(GraphGTest, testRandomEdge) {
	// we only test the uniform version
	count n = 4;
	count m = 5;
	count samples = 100000;
	double maxAbsoluteError = 0.005;

	Graph G = createGraph(n);
	G.addEdge(0, 1); // 0 * 1 = 0
	G.addEdge(1, 2); // 1 * 2 = 2
	G.addEdge(3, 2); // 3 * 2 = 1 (mod 5)
	G.addEdge(2, 2); // 2 * 2 = 4
	G.addEdge(3, 1); // 3 * 1 = 3
	ASSERT_EQ(m, G.numberOfEdges());

	std::vector<count> drawCounts(m, 0);
	for (auto e : G.randomEdges(samples)) {
		count id = (e.first * e.second) % 5;
		drawCounts[id]++;
	}
	for (node id = 0; id < m; id++) {
		double p = drawCounts[id] / (double)samples;
		ASSERT_NEAR(1.0 / m, p, maxAbsoluteError);
	}
}

/** GLOBAL PROPERTIES **/

TEST_P(GraphGTest, testSelfLoopCountSimple) {
	Graph G(Ghouse);
	G.addEdge(0, 0);
	EXPECT_EQ(1, G.numberOfSelfLoops());
}

TEST_P(GraphGTest, testIsWeighted) {
	ASSERT_EQ(isWeighted(), this->Ghouse.isWeighted());
}

TEST_P(GraphGTest, testIsDirected) {
	ASSERT_EQ(isDirected(), this->Ghouse.isDirected());
}

TEST_P(GraphGTest, testIsEmpty) {
	Graph G1 = createGraph(0);
	Graph G2 = createGraph(2);

	ASSERT_TRUE(G1.isEmpty());
	ASSERT_FALSE(G2.isEmpty());

	node v = G1.addNode();
	G2.removeNode(G2.randomNode());
	ASSERT_FALSE(G1.isEmpty());
	ASSERT_FALSE(G2.isEmpty());

	G1.removeNode(v);
	G2.removeNode(G2.randomNode());
	ASSERT_TRUE(G1.isEmpty());
	ASSERT_TRUE(G2.isEmpty());
}

TEST_P(GraphGTest, testNumberOfNodes) {
	ASSERT_EQ(this->n_house, this->Ghouse.numberOfNodes());

	Graph G1 = createGraph(0);
	ASSERT_EQ(0u, G1.numberOfNodes());
	G1.addNode();
	ASSERT_EQ(1u, G1.numberOfNodes());
	G1.addNode();
	ASSERT_EQ(2u, G1.numberOfNodes());
	G1.removeNode(0);
	ASSERT_EQ(1u, G1.numberOfNodes());
	G1.removeNode(1);
	ASSERT_EQ(0u, G1.numberOfNodes());
}

TEST_P(GraphGTest, testNumberOfEdges) {
	ASSERT_EQ(this->m_house, this->Ghouse.numberOfEdges());

	Graph G1 = createGraph(5);
	ASSERT_EQ(0u, G1.numberOfEdges());
	G1.addEdge(0, 1);
	ASSERT_EQ(1u, G1.numberOfEdges());
	G1.addEdge(1, 2);
	ASSERT_EQ(2u, G1.numberOfEdges());
	G1.removeEdge(0, 1);
	ASSERT_EQ(1u, G1.numberOfEdges());
	G1.removeEdge(1, 2);
	ASSERT_EQ(0u, G1.numberOfEdges());
}

TEST_P(GraphGTest, testNumberOfSelfLoops) {
	Graph G = createGraph(3);
	G.addEdge(0, 1);
	ASSERT_EQ(0u, G.numberOfSelfLoops());
	G.addEdge(0, 0);
	ASSERT_EQ(1u, G.numberOfSelfLoops());
	G.addEdge(1, 1);
	G.addEdge(1, 2);
	ASSERT_EQ(2u, G.numberOfSelfLoops());
	G.removeEdge(0, 0);
	ASSERT_EQ(1u, G.numberOfSelfLoops());
}

TEST_P(GraphGTest, testSelfLoopConversion) {
	Aux::Random::setSeed(1, false);
	const count runs = 100;
	const count n_max = 200;
	for (index i = 0; i < runs; i++) {
		bool directed = Aux::Random::probability() < 0.5;
		count n = Aux::Random::integer(n_max);
		Graph G(n, false, directed);

		G.forNodes([&](node v) {
			double p = Aux::Random::probability();

			if (p < 0.1) { // new node
				n++;
				G.addNode();
			} else {                                   // new edge
				node u = Aux::Random::integer(v, n - 1); // self-loops possible
				G.addEdge(v, u);
			}
		});
		count measuredSelfLoops = countSelfLoopsManually(G);
		EXPECT_EQ(G.numberOfSelfLoops(), measuredSelfLoops);
		Graph G_converted(G, false, !directed);
		EXPECT_EQ(G_converted.numberOfSelfLoops(), measuredSelfLoops);
	}
}

TEST_P(GraphGTest, testUpperNodeIdBound) {
	ASSERT_EQ(5u, this->Ghouse.upperNodeIdBound());

	Graph G1 = createGraph(0);
	ASSERT_EQ(0u, G1.upperNodeIdBound());
	G1.addNode();
	ASSERT_EQ(1u, G1.upperNodeIdBound());
	G1.addNode();
	ASSERT_EQ(2u, G1.upperNodeIdBound());
	G1.removeNode(1);
	ASSERT_EQ(2u, G1.upperNodeIdBound());
	G1.addNode();
	ASSERT_EQ(3u, G1.upperNodeIdBound());
}

TEST_P(GraphGTest, testCheckConsistency_MultiEdgeDetection) {
	Graph G = createGraph(3);
	ASSERT_TRUE(G.checkConsistency());
	G.addEdge(0, 1);
	ASSERT_TRUE(G.checkConsistency());
	G.addEdge(0, 2);
	G.addEdge(0, 1);
	ASSERT_FALSE(G.checkConsistency());
	G.removeEdge(0, 1);
	ASSERT_TRUE(G.checkConsistency());
	G.removeEdge(0, 1);
	ASSERT_TRUE(G.checkConsistency());
}

/** DYNAMICS **/

TEST_P(GraphGTest, testTime) {
	ASSERT_EQ(0u, this->Ghouse.time());
	this->Ghouse.timeStep();
	ASSERT_EQ(1u, this->Ghouse.time());
	this->Ghouse.timeStep();
	this->Ghouse.timeStep();
	ASSERT_EQ(3u, this->Ghouse.time());
}

/** EDGE ATTRIBUTES **/

TEST_P(GraphGTest, testWeight) {
	this->Ghouse.forNodes([&](node u) {
		this->Ghouse.forNodes([&](node v) {
			ASSERT_EQ(this->Ahouse[u][v], this->Ghouse.weight(u, v));
		});
	});
}

TEST_P(GraphGTest, testSetWeight) {
	Graph G = createGraph(10);
	G.addEdge(0, 1);
	G.addEdge(1, 2);

	if (isWeighted()) {
		// edges should get weight defaultWeight on creation and setWeight should
		// overwrite this
		G.setWeight(1, 2, 2.718);
		EXPECT_EQ(defaultEdgeWeight, G.weight(0, 1));
		EXPECT_EQ(2.718, G.weight(1, 2));
		if (isDirected()) {
			EXPECT_EQ(nullWeight, G.weight(1, 0));
			EXPECT_EQ(nullWeight, G.weight(2, 1));
		} else {
			// undirected graph is symmetric
			EXPECT_EQ(defaultEdgeWeight, G.weight(1, 0));
			EXPECT_EQ(2.718, G.weight(2, 1));
		}

		// setting an edge weight should create the edge if it doesn't exists
		ASSERT_FALSE(G.hasEdge(5, 6));
		G.setWeight(5, 6, 56.0);
		ASSERT_EQ(56.0, G.weight(5, 6));
		ASSERT_EQ(isDirected() ? nullWeight : 56.0, G.weight(6, 5));
		ASSERT_TRUE(G.hasEdge(5, 6));

		// directed graphs are not symmetric, undirected are
		G.setWeight(2, 1, 5.243);
		if (isDirected()) {
			EXPECT_EQ(2.718, G.weight(1, 2));
			EXPECT_EQ(5.243, G.weight(2, 1));
		} else {
			EXPECT_EQ(5.243, G.weight(1, 2));
			EXPECT_EQ(5.243, G.weight(2, 1));
		}

		// self-loop
		G.addEdge(4, 4, 2.5);
		ASSERT_EQ(2.5, G.weight(4, 4));
		G.setWeight(4, 4, 3.14);
		ASSERT_EQ(3.14, G.weight(4, 4));
	} else {
		EXPECT_ANY_THROW(G.setWeight(0, 1, 1.5));
	}
}

TEST_P(GraphGTest, increaseWeight) {
	Graph G = createGraph(5);
	G.addEdge(0, 1);
	G.addEdge(1, 2);
	G.addEdge(3, 4, 3.14);

	if (G.isWeighted()) {
		G.increaseWeight(1, 2, 0.5);
		G.increaseWeight(3, 4, -0.5);

		ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
		ASSERT_EQ(defaultEdgeWeight + 0.5, G.weight(1, 2));
		ASSERT_EQ(3.14 - 0.5, G.weight(3, 4));

		if (G.isDirected()) {
			// reverse edges do net exist => weight should be nullWeight
			ASSERT_EQ(nullWeight, G.weight(1, 0));
			ASSERT_EQ(nullWeight, G.weight(2, 1));
			ASSERT_EQ(nullWeight, G.weight(4, 3));
		} else {
			ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
			ASSERT_EQ(defaultEdgeWeight + 0.5, G.weight(2, 1));
			ASSERT_EQ(3.14 - 0.5, G.weight(3, 4));
		}
	} else {
		EXPECT_ANY_THROW(G.increaseWeight(1, 2, 0.3));
		EXPECT_ANY_THROW(G.increaseWeight(2, 3, 0.3)); // edge does not exists
	}
}

/** SUMS **/

TEST_P(GraphGTest, testTotalEdgeWeight) {
	Graph G1 = createGraph(5);
	Graph G2 = createGraph(5);
	G2.addEdge(0, 1, 3.14);

	if (this->Ghouse.isWeighted()) {
		ASSERT_EQ(0.0, G1.totalEdgeWeight());
		ASSERT_EQ(3.14, G2.totalEdgeWeight());
		ASSERT_EQ(36.0, this->Ghouse.totalEdgeWeight());
	} else {
		ASSERT_EQ(0 * defaultEdgeWeight, G1.totalEdgeWeight());
		ASSERT_EQ(1 * defaultEdgeWeight, G2.totalEdgeWeight());
		ASSERT_EQ(8 * defaultEdgeWeight, this->Ghouse.totalEdgeWeight());
	}
}

/** Collections **/

TEST_P(GraphGTest, testNodes) {
	Graph G = createGraph(3);
	G.addNode();
	G.removeNode(2);
	G.addNode();
	auto nodes = G.nodes();

	auto containsNode = [&nodes](node v) {
		return std::find(nodes.begin(), nodes.end(), v) != nodes.end();
	};

	ASSERT_EQ(G.numberOfNodes(), nodes.size());
	for (node v : nodes) {
		ASSERT_TRUE(containsNode(v));
	}
}

TEST_P(GraphGTest, testNeighbors) {
	auto neighbors = this->Ghouse.neighbors(1);
	auto containsNode = [&neighbors](node v) {
		return std::find(neighbors.begin(), neighbors.end(), v) != neighbors.end();
	};
	if (this->Ghouse.isDirected()) {
		ASSERT_TRUE(containsNode(0));
		ASSERT_TRUE(containsNode(4));
	} else {
		ASSERT_TRUE(containsNode(0));
		ASSERT_TRUE(containsNode(2));
		ASSERT_TRUE(containsNode(3));
		ASSERT_TRUE(containsNode(4));
	}
}

TEST_P(GraphGTest, testEdges) {
	// add self-loop
	this->Ghouse.addEdge(3, 3);
	auto isCorrectEdge = [&](node u, node v) {
		if (u == 3 && v == 3) {
			return true;
		}
		auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
		                    std::make_pair(u, v));
		if (it != this->houseEdgesOut.end()) {
			return true;
		} else if (!this->Ghouse.isDirected()) {
			it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
			               std::make_pair(v, u));
			return it != this->houseEdgesOut.end();
		}
		return false;
	};

	auto edges = this->Ghouse.edges();
	ASSERT_EQ(this->m_house + 1, edges.size()); // plus self-loop
	for (auto e : edges) {
		ASSERT_TRUE(isCorrectEdge(e.first, e.second))
		    << "(" << e.first << ", " << e.second
		    << ") is in edge array, but is not an edge of Ghouse";
	}
}

TEST_P(GraphGTest, testTranspose) {
	Graph G = this->Ghouse;

	G.addNode(); // node 5
	G.addNode(); // node 6
	G.removeNode(5);

	if (!G.isWeighted()) {
		G.addEdge(0, 0);
		G.addEdge(0, 4);
		G.removeEdge(0, 4);
		G.addEdge(0, 6);
	} else {
		G.addEdge(0, 0, 3.14);
		G.addEdge(0, 4, 3.14);
		G.removeEdge(0, 4);
		G.addEdge(0, 6, 3.14);
	}
	DEBUG("Ghouse: ", G.nodes(), " ", G.edges());
	// expect throw error when G is undirected
	if (!G.isDirected()) {
		EXPECT_ANY_THROW(G.transpose());
	} else {
		Graph Gtrans = G.transpose();
		DEBUG("GhouseTrans: ", Gtrans.nodes(), " ", Gtrans.edges());
		// check summation statistics
		EXPECT_EQ(G.numberOfNodes(), Gtrans.numberOfNodes());
		EXPECT_EQ(G.numberOfEdges(), Gtrans.numberOfEdges());
		EXPECT_EQ(G.totalEdgeWeight(), Gtrans.totalEdgeWeight());
		EXPECT_EQ(G.numberOfSelfLoops(), Gtrans.numberOfSelfLoops());

		// check time step
		EXPECT_EQ(G.time(), Gtrans.time());

		// check graph names
		EXPECT_EQ(G.getName() + "Transpose", Gtrans.getName());

		// test for regular edges
		EXPECT_TRUE(G.hasEdge(0, 6));
		EXPECT_FALSE(G.hasEdge(6, 0));
		EXPECT_TRUE(Gtrans.hasEdge(6, 0));
		EXPECT_FALSE(Gtrans.hasEdge(0, 6));
		// .. and for selfloops
		EXPECT_TRUE(G.hasEdge(0, 0));
		EXPECT_TRUE(Gtrans.hasEdge(0, 0));

		// check for edge weights
		if (G.isWeighted()) {
			EXPECT_EQ(G.weight(0, 6), 3.14);
			EXPECT_EQ(Gtrans.weight(6, 0), 3.14);
			EXPECT_EQ(G.weight(0, 0), Gtrans.weight(0, 0));
		}
	}
}

/** NODE ITERATORS **/

TEST_P(GraphGTest, testForNodes) {
	Graph G = createGraph(3);
	std::vector<bool> visited(4, false);
	G.forNodes([&](node v) {
		ASSERT_FALSE(visited[v]);
		if (v == 2) {
			G.addNode();
		}
		visited[v] = true;
	});
	for (bool b : visited) {
		ASSERT_TRUE(b);
	}
}

TEST_P(GraphGTest, testParallelForNodes) {
	std::vector<node> visited;
	this->Ghouse.parallelForNodes([&](node u) {
#pragma omp critical
		visited.push_back(u);
	});

	Aux::Parallel::sort(visited.begin(), visited.end());

	ASSERT_EQ(5u, visited.size());
	for (index i = 0; i < this->Ghouse.upperNodeIdBound(); i++) {
		ASSERT_EQ(i, visited[i]);
	}
}

TEST_P(GraphGTest, forNodesWhile) {
	count n = 100;
	Graph G = createGraph(n);
	count stopAfter = 10;
	count nodesSeen = 0;

	G.forNodesWhile([&]() { return nodesSeen < stopAfter; },
	                [&](node u) { nodesSeen++; });

	ASSERT_EQ(stopAfter, nodesSeen);
}

TEST_P(GraphGTest, testForNodesInRandomOrder) {
	count n = 1000;
	count samples = 100;
	double maxAbsoluteError = 0.005;
	Graph G = createGraph(n);

	node lastNode = n / 2;
	count greaterLastNode = 0;
	count smallerLastNode = 0;
	std::vector<count> visitCount(n, 0);

	for (count i = 0; i < samples; i++) {
		G.forNodesInRandomOrder([&](node v) {
			if (v > lastNode) {
				greaterLastNode++;
			} else {
				smallerLastNode++;
			}
			visitCount[v]++;
			lastNode = v;
		});
	}

	for (node v = 0; v < n; v++) {
		ASSERT_EQ(samples, visitCount[v]);
	}

	ASSERT_NEAR(0.5, (double)greaterLastNode / n / samples, maxAbsoluteError);
	ASSERT_NEAR(0.5, (double)smallerLastNode / n / samples, maxAbsoluteError);
}

TEST_P(GraphGTest, testForNodePairs) {
	count n = 10;
	count m = n * (n - 1) / 2;
	Graph G = createGraph(n);

	// add all edges
	G.forNodePairs([&](node u, node v) {
		ASSERT_FALSE(G.hasEdge(u, v));
		G.addEdge(u, v);
		ASSERT_TRUE(G.hasEdge(u, v));
	});

	EXPECT_EQ(m, G.numberOfEdges());

	// remove all edges
	G.forNodePairs([&](node u, node v) {
		ASSERT_TRUE(G.hasEdge(u, v));
		G.removeEdge(u, v);
		ASSERT_FALSE(G.hasEdge(u, v));
	});

	EXPECT_EQ(0u, G.numberOfEdges());
}

/** EDGE ITERATORS **/

TEST_P(GraphGTest, testForEdges) {
	Graph G = createGraph(4);
	G.addEdge(0, 1); // 0 * 1 = 0
	G.addEdge(1, 2); // 1 * 2 = 2
	G.addEdge(3, 2); // 3 * 2 = 1 (mod 5)
	G.addEdge(2, 2); // 2 * 2 = 4
	G.addEdge(3, 1); // 3 * 1 = 3

	std::vector<bool> edgesSeen(5, false);

	G.forEdges([&](node u, node v) {
		ASSERT_TRUE(G.hasEdge(u, v));
		index id = (u * v) % 5;
		edgesSeen[id] = true;
	});

	for (auto b : edgesSeen) {
		ASSERT_TRUE(b);
	}
}

TEST_P(GraphGTest, testForWeightedEdges) {
	double epsilon = 1e-6;

	Graph G = createGraph(4);
	G.addEdge(0, 1, 0.1); // 0 * 1 = 0
	G.addEdge(3, 2, 0.2); // 3 * 2 = 1 (mod 5)
	G.addEdge(1, 2, 0.3); // 1 * 2 = 2
	G.addEdge(3, 1, 0.4); // 3 * 1 = 3
	G.addEdge(2, 2, 0.5); // 2 * 2 = 4

	std::vector<bool> edgesSeen(5, false);

	edgeweight weightSum = 0;
	G.forEdges([&](node u, node v, edgeweight ew) {
		ASSERT_TRUE(G.hasEdge(u, v));
		ASSERT_EQ(G.weight(u, v), ew);

		index id = (u * v) % 5;
		edgesSeen[id] = true;
		if (G.isWeighted()) {
			ASSERT_NEAR((id + 1) * 0.1, ew, epsilon);
		} else {
			ASSERT_EQ(defaultEdgeWeight, ew);
		}
		weightSum += ew;
	});

	for (auto b : edgesSeen) {
		ASSERT_TRUE(b);
	}
	if (G.isWeighted()) {
		ASSERT_NEAR(1.5, weightSum, epsilon);
	} else {
		ASSERT_NEAR(5 * defaultEdgeWeight, weightSum, epsilon);
	}
}

TEST_P(GraphGTest, testParallelForWeightedEdges) {
	count n = 4;
	Graph G = createGraph(n);
	G.forNodePairs([&](node u, node v) { G.addEdge(u, v, 1.0); });

	edgeweight weightSum = 0.0;
	G.parallelForEdges([&](node u, node v, edgeweight ew) {
#pragma omp atomic
		weightSum += ew;
	});

	ASSERT_EQ(6.0, weightSum) << "sum of edge weights should be 6 in every case";
}

TEST_P(GraphGTest, testParallelForEdges) {
	count n = 4;
	Graph G = createGraph(n);
	G.forNodePairs([&](node u, node v) { G.addEdge(u, v); });

	edgeweight weightSum = 0.0;
	G.parallelForEdges([&](node u, node v) {
#pragma omp atomic
		weightSum += 1;
	});

	ASSERT_EQ(6.0, weightSum) << "sum of edge weights should be 6 in every case";
}

// template<typename L> void forEdgesWithAttribute_double(int attrId, L handle)
// const;

/** NEIGHBORHOOD ITERATORS **/

TEST_P(GraphGTest, testForNeighborsOf) {
	std::vector<node> visited;
	this->Ghouse.forNeighborsOf(3, [&](node u) { visited.push_back(u); });

	Aux::Parallel::sort(visited.begin(), visited.end());

	if (isDirected()) {
		ASSERT_EQ(2u, visited.size());
		ASSERT_EQ(1u, visited[0]);
		ASSERT_EQ(2u, visited[1]);
	} else {
		ASSERT_EQ(3u, visited.size());
		ASSERT_EQ(1u, visited[0]);
		ASSERT_EQ(2u, visited[1]);
		ASSERT_EQ(4u, visited[2]);
	}
}

TEST_P(GraphGTest, testForWeightedNeighborsOf) {
	std::vector<std::pair<node, edgeweight>> visited;
	this->Ghouse.forNeighborsOf(3, [&](node u, edgeweight ew) {
		visited.push_back(std::make_pair(u, ew));
	});

	// should sort after the first element
	Aux::Parallel::sort(visited.begin(), visited.end());

	if (isGraph()) {
		ASSERT_EQ(3u, visited.size());
		ASSERT_EQ(1u, visited[0].first);
		ASSERT_EQ(2u, visited[1].first);
		ASSERT_EQ(4u, visited[2].first);
		ASSERT_EQ(defaultEdgeWeight, visited[0].second);
		ASSERT_EQ(defaultEdgeWeight, visited[1].second);
		ASSERT_EQ(defaultEdgeWeight, visited[2].second);
	}

	if (isWeightedGraph()) {
		ASSERT_EQ(3u, visited.size());
		ASSERT_EQ(1u, visited[0].first);
		ASSERT_EQ(2u, visited[1].first);
		ASSERT_EQ(4u, visited[2].first);
		ASSERT_EQ(1.0, visited[0].second);
		ASSERT_EQ(7.0, visited[1].second);
		ASSERT_EQ(6.0, visited[2].second);
	}

	if (isDirectedGraph()) {
		ASSERT_EQ(2u, visited.size());
		ASSERT_EQ(1u, visited[0].first);
		ASSERT_EQ(2u, visited[1].first);
		ASSERT_EQ(defaultEdgeWeight, visited[0].second);
		ASSERT_EQ(defaultEdgeWeight, visited[1].second);
	}

	if (isWeightedDirectedGraph()) {
		ASSERT_EQ(2u, visited.size());
		ASSERT_EQ(1u, visited[0].first);
		ASSERT_EQ(2u, visited[1].first);
		ASSERT_EQ(1.0, visited[0].second);
		ASSERT_EQ(7.0, visited[1].second);
	}
}

TEST_P(GraphGTest, testForEdgesOf) {
	count m = 0;
	std::vector<int> visited(this->m_house, 0);

	this->Ghouse.forNodes([&](node u) {
		this->Ghouse.forEdgesOf(u, [&](node v, node w) {
			// edges should be v to w, so if we iterate over edges from u, u should be
			// equal v
			EXPECT_EQ(u, v);

			auto e = std::make_pair(v, w);
			auto it =
			    std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
			if (!isDirected() && it == this->houseEdgesOut.end()) {
				auto e2 = std::make_pair(w, v);
				it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
				               e2);
			}

			EXPECT_TRUE(it != this->houseEdgesOut.end());

			// find index in edge array
			int i = std::distance(this->houseEdgesOut.begin(), it);
			if (isDirected()) {
				// make sure edge was not visited before (would be visited twice)
				EXPECT_EQ(0, visited[i]);
			}

			// mark edge as visited
			visited[i]++;
			m++;
		});
	});

	if (isDirected()) {
		// we iterated over all outgoing edges once
		EXPECT_EQ(this->m_house, m);
		for (auto c : visited) {
			EXPECT_EQ(1, c);
		}
	} else {
		// we iterated over all edges in both directions
		EXPECT_EQ(2 * this->m_house, m);
		for (auto c : visited) {
			EXPECT_EQ(2, c);
		}
	}
}

TEST_P(GraphGTest, testForWeightedEdgesOf) {
	count m = 0;
	std::vector<int> visited(this->m_house, 0);
	double sumOfWeights = 0;

	this->Ghouse.forNodes([&](node u) {
		this->Ghouse.forEdgesOf(u, [&](node v, node w, edgeweight ew) {
			// edges should be v to w, so if we iterate over edges from u, u should be
			// equal v
			EXPECT_EQ(u, v);
			sumOfWeights += ew;
			auto e = std::make_pair(v, w);
			auto it =
			    std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
			if (!isDirected() && it == this->houseEdgesOut.end()) {
				auto e2 = std::make_pair(w, v);
				it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(),
				               e2);
			}

			EXPECT_TRUE(it != this->houseEdgesOut.end());

			// find index in edge array
			int i = std::distance(this->houseEdgesOut.begin(), it);
			if (isDirected()) {
				// make sure edge was not visited before (would be visited twice)
				EXPECT_EQ(0, visited[i]);
			}

			// mark edge as visited
			visited[i]++;
			m++;
		});
	});

	if (isGraph()) {
		EXPECT_EQ(sumOfWeights, m);
		EXPECT_EQ(2 * this->m_house, m);
		for (auto c : visited) {
			EXPECT_EQ(2, c);
		}
	}

	if (isWeightedGraph()) {
		// we iterated over all edges in both directions
		EXPECT_EQ(2 * this->m_house, m);
		EXPECT_EQ(sumOfWeights, 72);
		for (auto c : visited) {
			EXPECT_EQ(2, c);
		}
	}

	if (isDirectedGraph()) {
		// we iterated over all outgoing edges once
		EXPECT_EQ(this->m_house, m);
		EXPECT_EQ(sumOfWeights, m);
		for (auto c : visited) {
			EXPECT_EQ(1, c);
		}
	}

	if (isWeightedDirectedGraph()) {
		EXPECT_EQ(sumOfWeights, 36);
		EXPECT_EQ(this->m_house, m);
		for (auto c : visited) {
			EXPECT_EQ(1, c);
		}
	}
}

TEST_P(GraphGTest, testForInNeighborsOf) {
	std::vector<node> visited;
	this->Ghouse.forInNeighborsOf(2, [&](node v) { visited.push_back(v); });
	Aux::Parallel::sort(visited.begin(), visited.end());

	if (isDirected()) {
		EXPECT_EQ(2u, visited.size());
		EXPECT_EQ(0u, visited[0]);
		EXPECT_EQ(3u, visited[1]);
	} else {
		EXPECT_EQ(4u, visited.size());
		EXPECT_EQ(0u, visited[0]);
		EXPECT_EQ(1u, visited[1]);
		EXPECT_EQ(3u, visited[2]);
		EXPECT_EQ(4u, visited[3]);
	}
}

TEST_P(GraphGTest, testForWeightedInNeighborsOf) {
	std::vector<std::pair<node, edgeweight>> visited;
	this->Ghouse.forInNeighborsOf(3, [&](node v, edgeweight ew) {
		visited.push_back({v, ew});
	});
	Aux::Parallel::sort(visited.begin(), visited.end());

	if (isGraph()) {
		ASSERT_EQ(3u, visited.size());
		ASSERT_EQ(1u, visited[0].first);
		ASSERT_EQ(2u, visited[1].first);
		ASSERT_EQ(4u, visited[2].first);
		ASSERT_EQ(defaultEdgeWeight, visited[0].second);
		ASSERT_EQ(defaultEdgeWeight, visited[1].second);
		ASSERT_EQ(defaultEdgeWeight, visited[2].second);
	}

	if (isWeightedGraph()) {
		ASSERT_EQ(3u, visited.size());
		ASSERT_EQ(1u, visited[0].first);
		ASSERT_EQ(2u, visited[1].first);
		ASSERT_EQ(4u, visited[2].first);
		ASSERT_EQ(1.0, visited[0].second);
		ASSERT_EQ(7.0, visited[1].second);
		ASSERT_EQ(6.0, visited[2].second);
	}

	if (isDirectedGraph()) {
		ASSERT_EQ(1u, visited.size());
		ASSERT_EQ(4u, visited[0].first);
		ASSERT_EQ(defaultEdgeWeight, visited[0].second);
	}

	if (isWeightedDirectedGraph()) {
		ASSERT_EQ(1u, visited.size());
		ASSERT_EQ(4u, visited[0].first);
		ASSERT_EQ(6.0, visited[0].second);
	}
}

TEST_P(GraphGTest, forInEdgesOf) {
	std::vector<bool> visited(this->n_house, false);
	this->Ghouse.forInEdgesOf(3, [&](node u, node v) {
		ASSERT_EQ(3u, v);
		ASSERT_TRUE(this->Ahouse[u][v] > 0.0);
		ASSERT_TRUE(this->Ghouse.hasEdge(u, v));
		ASSERT_FALSE(visited[u]);
		visited[u] = true;
	});

	if (isDirected()) {
		EXPECT_FALSE(visited[0]);
		EXPECT_FALSE(visited[1]);
		EXPECT_FALSE(visited[2]);
		EXPECT_FALSE(visited[3]);
		EXPECT_TRUE(visited[4]);
	} else {
		EXPECT_FALSE(visited[0]);
		EXPECT_TRUE(visited[1]);
		EXPECT_TRUE(visited[2]);
		EXPECT_FALSE(visited[3]);
		EXPECT_TRUE(visited[4]);
	}
}

TEST_P(GraphGTest, testForWeightedInEdgesOf) {
	// add self-loop
	this->Ghouse.addEdge(3, 3, 2.5);
	this->Ahouse[3][3] = 2.5;

	std::vector<edgeweight> visited(this->n_house, -1.0);
	this->Ghouse.forInEdgesOf(3, [&](node v, node u, edgeweight ew) {
		ASSERT_EQ(-1.0, visited[v]);
		visited[u] = ew;
	});

	if (isGraph()) {
		ASSERT_EQ(-1.0, visited[0]);
		ASSERT_EQ(defaultEdgeWeight, visited[1]);
		ASSERT_EQ(defaultEdgeWeight, visited[2]);
		ASSERT_EQ(defaultEdgeWeight, visited[3]);
		ASSERT_EQ(defaultEdgeWeight, visited[4]);
	}

	if (isWeightedGraph()) {
		ASSERT_EQ(-1.0, visited[0]);
		ASSERT_EQ(this->Ahouse[3][1], visited[1]);
		ASSERT_EQ(this->Ahouse[3][2], visited[2]);
		ASSERT_EQ(this->Ahouse[3][3], visited[3]);
		ASSERT_EQ(this->Ahouse[3][4], visited[4]);
	}

	if (isDirectedGraph()) {
		ASSERT_EQ(-1.0, visited[0]);
		ASSERT_EQ(-1.0, visited[1]);
		ASSERT_EQ(-1.0, visited[2]);
		ASSERT_EQ(defaultEdgeWeight, visited[3]);
		ASSERT_EQ(defaultEdgeWeight, visited[4]);
	}

	if (isWeightedDirectedGraph()) {
		ASSERT_EQ(-1.0, visited[0]);
		ASSERT_EQ(-1.0, visited[1]);
		ASSERT_EQ(-1.0, visited[2]);
		ASSERT_EQ(this->Ahouse[3][3], visited[3]);
		ASSERT_EQ(this->Ahouse[4][3], visited[4]);
	}
}

/** REDUCTION ITERATORS **/

TEST_P(GraphGTest, testParallelSumForNodes) {
	count n = 10;
	Graph G = createGraph(n);
	double sum = G.parallelSumForNodes([](node v) { return 2 * v + 0.5; });

	double expected_sum = n * (n - 1) + n * 0.5;
	ASSERT_EQ(expected_sum, sum);
}

TEST_P(GraphGTest, testParallelSumForWeightedEdges) {
	double sum = this->Ghouse.parallelSumForEdges(
	    [](node u, node v, edgeweight ew) { return 1.5 * ew; });

	double expected_sum = 1.5 * this->Ghouse.totalEdgeWeight();
	ASSERT_EQ(expected_sum, sum);
}

/** GRAPH SEARCHES **/

TEST_P(GraphGTest, testBFSfrom) {
	std::vector<count> visitedOrder(5, none);
	index i = 0;
	this->Ghouse.BFSfrom(3, [&](node v, count dist) {
		EXPECT_EQ(none, visitedOrder[v]); // visit every node once
		visitedOrder[v] = i++;
	});
	// have we visited all nodes
	for (count l : visitedOrder) {
		EXPECT_TRUE(l != none);
	}

	if (isDirected()) {
		// root on level 0
		EXPECT_EQ(0u, visitedOrder[3]);

		// level 1
		EXPECT_TRUE((visitedOrder[1] == 1) ^ (visitedOrder[1] == 2));
		EXPECT_TRUE((visitedOrder[2] == 1) ^ (visitedOrder[2] == 2));

		// level 2
		EXPECT_TRUE((visitedOrder[0] == 3) ^ (visitedOrder[0] == 4));
		EXPECT_TRUE((visitedOrder[4] == 3) ^ (visitedOrder[4] == 4));
	} else {
		EXPECT_EQ(0u, visitedOrder[3]);
		EXPECT_TRUE((visitedOrder[1] == 1) ^ (visitedOrder[1] == 2) ^
		            (visitedOrder[1] == 3));
		EXPECT_TRUE((visitedOrder[2] == 1) ^ (visitedOrder[2] == 2) ^
		            (visitedOrder[2] == 3));
		EXPECT_TRUE((visitedOrder[4] == 1) ^ (visitedOrder[4] == 2) ^
		            (visitedOrder[4] == 3));
		EXPECT_TRUE((visitedOrder[0] == 4));
	}
}

TEST_P(GraphGTest, testDFSfrom) {
	if (isDirected()) {
		std::vector<count> visitedOrder(5, none);
		index i = 0;
		this->Ghouse.DFSfrom(3, [&](node v) {
			EXPECT_EQ(none, visitedOrder[v]); // visit every node once
			visitedOrder[v] = i++;
		});

		// have we visited all nodes
		for (count l : visitedOrder) {
			EXPECT_TRUE(l != none);
		}

		// root on level 0
		EXPECT_EQ(0u, visitedOrder[3]);

		// level 1
		EXPECT_TRUE((visitedOrder[1] == 1) ^ (visitedOrder[2] == 1));

		// level 2
		EXPECT_TRUE((visitedOrder[0] == 2) ^ (visitedOrder[1] == 2) ^
		            (visitedOrder[4] == 2));

		// level 3
		EXPECT_TRUE((visitedOrder[2] == 3) ^ (visitedOrder[0] == 3) ^
		            (visitedOrder[4] == 3) ^ (visitedOrder[1] == 3));

		// level 4
		EXPECT_TRUE((visitedOrder[2] == 4) ^ (visitedOrder[4] == 4) ^
		            (visitedOrder[0] == 4));
	} else {
		count n = 5;
		std::vector<count> visitedOrder(n, none);
		Graph G = createGraph(n);
		G.addEdge(0, 1);
		G.addEdge(0, 2);
		G.addEdge(2, 3);
		G.addEdge(3, 4);

		index i = 0;
		G.DFSfrom(0, [&](node v) { visitedOrder[v] = i++; });

		for (count l : visitedOrder) {
			EXPECT_TRUE(l != none);
		}

		EXPECT_EQ(0u, visitedOrder[0]);

		EXPECT_TRUE((visitedOrder[1] == 1) ^ (visitedOrder[2] == 1));

		EXPECT_TRUE((visitedOrder[2] == 2) ^ (visitedOrder[3] == 2));

		EXPECT_TRUE((visitedOrder[3] == 3) ^ (visitedOrder[4] == 3));

		EXPECT_TRUE((visitedOrder[4] == 4) ^ (visitedOrder[1] == 4));
	}
}

TEST_P(GraphGTest, testEdgeIndexGenerationDirected) {
	Graph G = Graph(10, false, true);
	G.addEdge(2, 0);
	G.addEdge(2, 1);
	G.addEdge(2, 2);
	G.addEdge(5, 6);
	G.addEdge(6, 5);
	G.addEdge(1, 2);

	G.indexEdges();

	// Check consecutiveness of edgeids according to edge iterators
	edgeid expectedId = 0;
	G.forEdges([&](node u, node v) { EXPECT_EQ(expectedId++, G.edgeId(u, v)); });

	// Add some more edges
	EXPECT_EQ(6, G.upperEdgeIdBound());
	G.addEdge(8, 9);
	EXPECT_EQ(7, G.upperEdgeIdBound());
	G.addEdge(9, 8);

	// Check that asymmetric edges do not have the same id
	EXPECT_NE(G.edgeId(6, 5), G.edgeId(5, 6));
	EXPECT_NE(G.edgeId(2, 1), G.edgeId(1, 2));
	EXPECT_NE(G.edgeId(9, 8), G.edgeId(8, 9));
	EXPECT_EQ(7, G.edgeId(9, 8));
	EXPECT_EQ(8, G.upperEdgeIdBound());
}

TEST_P(GraphGTest, testEdgeIndexGenerationUndirected) {
	Graph G = Graph(10, false, false);

	G.addEdge(0, 0);
	G.addEdge(2, 0);
	G.addEdge(2, 1);
	G.addEdge(2, 2);
	G.addEdge(5, 6);

	G.indexEdges();

	// Check consecutiveness of edgeids according to edge iterators
	edgeid expectedId = 0;
	G.forEdges([&](node u, node v) { EXPECT_EQ(expectedId++, G.edgeId(u, v)); });

	// Add some more edges. This will likely destroy consecutiveness...
	G.addEdge(3, 4);
	G.addEdge(7, 8);
	EXPECT_EQ(6, G.edgeId(7, 8));
	EXPECT_EQ(7, G.upperEdgeIdBound());

	// Anyway, check uniqueness and validity of the edgeids
	std::set<edgeid> ids;
	edgeid upperEdgeIdBound = G.upperEdgeIdBound();

	G.forEdges([&](node u, node v) {
		edgeid id = G.edgeId(u, v);
		EXPECT_EQ(id, G.edgeId(v, u));
		EXPECT_LT(id, upperEdgeIdBound);

		EXPECT_NE(none, id);
		EXPECT_FALSE(ids.erase(id));
		ids.insert(id);
	});
}

TEST_P(GraphGTest, testForEdgesWithIds) {
	std::vector<Graph> graphs;
	graphs.push_back(Graph(10, false, false));
	graphs.push_back(Graph(10, false, true));
	graphs.push_back(Graph(10, true, false));
	graphs.push_back(Graph(10, true, true));

	for (auto graph = graphs.begin(); graph != graphs.end(); ++graph) {
		graph->addEdge(0, 0);
		graph->addEdge(1, 2);
		graph->addEdge(4, 5);

		// No edge indices

		count m = 0;
		graph->forEdges([&](node u, node v, edgeid eid) {
			EXPECT_EQ(0, eid);
			m++;
		});
		ASSERT_EQ(3u, m);

		m = 0;
		graph->parallelForEdges([&](node u, node v, edgeid eid) {
			EXPECT_EQ(0, eid);
#pragma omp atomic
			m++;
		});
		ASSERT_EQ(3u, m);

		// With edge indices
		graph->indexEdges();

		edgeid expectedId = 0;
		m = 0;
		graph->forEdges([&](node u, node v, edgeid eid) {
			EXPECT_EQ(expectedId++, eid);
			EXPECT_LT(eid, graph->upperEdgeIdBound());
			m++;
		});
		ASSERT_EQ(3u, m);

		m = 0;
		graph->parallelForEdges([&](node u, node v, edgeid eid) {
			EXPECT_NE(none, eid);
			EXPECT_LT(eid, graph->upperEdgeIdBound());
#pragma omp atomic
			m++;
		});
		ASSERT_EQ(3u, m);
	}
}

TEST_P(GraphGTest, testForWeightedEdgesWithIds) {
	std::vector<Graph> graphs;
	graphs.push_back(Graph(10, false, false));
	graphs.push_back(Graph(10, false, true));
	graphs.push_back(Graph(10, true, false));
	graphs.push_back(Graph(10, true, true));

	for (auto graph = graphs.begin(); graph != graphs.end(); ++graph) {
		graph->addEdge(0, 0, 2);
		graph->addEdge(1, 2, 2);
		graph->addEdge(4, 5, 2);

		// No edge indices

		count m = 0;
		edgeweight sum = 0;
		graph->forEdges([&](node u, node v, edgeweight ew, edgeid eid) {
			EXPECT_EQ(0, eid);
			m++;
			sum += ew;
		});
		ASSERT_EQ(3u, m);

		if (graph->isWeighted()) {
			ASSERT_EQ(6.0, sum);
		} else {
			ASSERT_EQ(3.0, sum);
		}

		m = 0;
		sum = .0;
		graph->parallelForEdges([&](node u, node v, edgeweight ew, edgeid eid) {
			EXPECT_EQ(0, eid);
#pragma omp atomic
			m++;
#pragma omp atomic
			sum += ew;
		});
		ASSERT_EQ(3u, m);

		if (graph->isWeighted()) {
			ASSERT_EQ(6.0, sum);
		} else {
			ASSERT_EQ(3.0, sum);
		}

		// With edge indices
		graph->indexEdges();

		edgeid expectedId = 0;
		m = 0;
		sum = .0;
		graph->forEdges([&](node u, node v, edgeweight ew, edgeid eid) {
			EXPECT_EQ(expectedId++, eid);
			EXPECT_LT(eid, graph->upperEdgeIdBound());
			m++;
			sum += ew;
		});
		ASSERT_EQ(3u, m);

		if (graph->isWeighted()) {
			ASSERT_EQ(6.0, sum);
		} else {
			ASSERT_EQ(3.0, sum);
		}

		m = 0;
		sum = .0;
		graph->parallelForEdges([&](node u, node v, edgeweight ew, edgeid eid) {
			EXPECT_NE(none, eid);
			EXPECT_LT(eid, graph->upperEdgeIdBound());
#pragma omp atomic
			m++;
#pragma omp atomic
			sum += ew;
		});
		ASSERT_EQ(3u, m);

		if (graph->isWeighted()) {
			ASSERT_EQ(6.0, sum);
		} else {
			ASSERT_EQ(3.0, sum);
		}
	}
}

/*TEST_P(GraphGTest, testInForEdgesUndirected) {
  METISGraphReader reader;
  Graph G = reader.read("input/PGPgiantcompo.graph");
  DEBUG(G.upperNodeIdBound());
  node u = 5474;
  G.forInEdgesOf(u, [&](node u, node z, edgeweight w){
    DEBUG("(1) node: ", u, " neigh:", z, " weight: ", w);
  });
  G.forEdgesOf(u, [&](node u, node z, edgeweight w){
    DEBUG("(2) node: ", u, " neigh:", z, " weight: ", w);
  });


  node source = 1492;
  DynBFS bfs(G, source, false);
  bfs.run();

  std::vector<std::pair<node, double> > choices1;
  G.forInEdgesOf(5474, [&](node t, node z, edgeweight w){
    INFO("considered edge (1): ", t, z, w);
    if (Aux::NumericTools::logically_equal(bfs.distance(t), bfs.distance(z) +
w)) {
      // workaround for integer overflow in large graphs
      bigfloat tmp = bfs.getNumberOfPaths(z) / bfs.getNumberOfPaths(t);
      double weight;
      tmp.ToDouble(weight);
      choices1.emplace_back(z, weight);
    }
  });
  std::vector<std::pair<node, double> > choices2;
  G.forEdgesOf(5474, [&](node t, node z, edgeweight w){
    INFO("considered edge (2): ", t, z, w);
    if (Aux::NumericTools::logically_equal(bfs.distance(t), bfs.distance(z) +
w)) {
      // workaround for integer overflow in large graphs
      bigfloat tmp = bfs.getNumberOfPaths(z) / bfs.getNumberOfPaths(t);
      double weight;
      tmp.ToDouble(weight);
      choices2.emplace_back(z, weight);
    }
  });

  INFO(choices1);
  INFO(choices2);
}
*/

TEST_P(GraphGTest, testCompactEdges) {
	Graph G = this->Ghouse;
	G.indexEdges();

	G.addEdge(0, 4);
	G.addEdge(0, 3);
	G.removeEdge(0, 3);
	G.removeEdge(0, 4);

	G.compactEdges();

	std::vector<std::pair<node, node>> outEdges;
	outEdges.reserve(this->Ghouse.numberOfEdges());

	this->Ghouse.forEdges([&](node u, node v) { outEdges.emplace_back(u, v); });

	auto it = outEdges.begin();

	G.forEdges([&](node u, node v) {
		ASSERT_NE(it, outEdges.end());
		EXPECT_EQ(it->first, u);
		EXPECT_EQ(it->second, v);
		++it;
	});
}

TEST_P(GraphGTest, testSortEdges) {
	Graph G = this->Ghouse;

	for (int i = 0; i < 2; ++i) {
		if (i > 0) {
			G.indexEdges();
		}

		Graph origG = G;

		G.sortEdges();

		std::vector<std::tuple<node, node, edgeweight, edgeid>> edges;
		edges.reserve(origG.numberOfEdges() * 4);

		std::vector<std::tuple<node, edgeweight, edgeid>> outEdges;
		origG.forNodes([&](node u) {
			origG.forEdgesOf(u, [&](node u, node v, edgeweight w, edgeid eid) {
				outEdges.emplace_back(v, w, eid);
			});

			Aux::Parallel::sort(outEdges.begin(), outEdges.end());

			for (auto x : outEdges) {
				edges.emplace_back(u, std::get<0>(x), std::get<1>(x), std::get<2>(x));
			}

			outEdges.clear();

			origG.forInEdgesOf(u, [&](node u, node v, edgeweight w, edgeid eid) {
				outEdges.emplace_back(v, w, eid);
			});

			Aux::Parallel::sort(outEdges.begin(), outEdges.end());

			for (auto x : outEdges) {
				edges.emplace_back(u, std::get<0>(x), std::get<1>(x), std::get<2>(x));
			}

			outEdges.clear();
		});

		auto it = edges.begin();

		G.forNodes([&](node u) {
			G.forEdgesOf(u, [&](node u, node v, edgeweight w, edgeid eid) {
				ASSERT_NE(it, edges.end());
				EXPECT_EQ(*it, std::make_tuple(u, v, w, eid))
				    << "Out edge (" << u << ", " << v << ", " << w << ", " << eid
				    << ") was expected to be (" << std::get<0>(*it) << ", "
				    << std::get<1>(*it) << ", " << std::get<2>(*it) << ", "
				    << std::get<3>(*it) << ")";
				++it;
			});
			G.forInEdgesOf(u, [&](node u, node v, edgeweight w, edgeid eid) {
				ASSERT_NE(it, edges.end());
				EXPECT_EQ(*it, std::make_tuple(u, v, w, eid))
				    << "In edge (" << u << ", " << v << ", " << w << ", " << eid
				    << ") was expected to be (" << std::get<0>(*it) << ", "
				    << std::get<1>(*it) << ", " << std::get<2>(*it) << ", "
				    << std::get<3>(*it) << ")";
				++it;
			});
		});
	}
}

} /* namespace NetworKit */
