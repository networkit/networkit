/*
 * ConnectedComponentsGTest.cpp
 *
 *  Created on: Sep 16, 2013
 *      Author: Maximilian Vogel
 */
#ifndef NOGTEST

#include "ConnectedComponentsGTest.h"
#include "../ConnectedComponents.h"
#include "../ParallelConnectedComponents.h"
#include "../StronglyConnectedComponents.h"
#include "../GraphProperties.h"
#include "../Diameter.h"
#include "../../io/METISGraphReader.h"
#include "../../generators/HavelHakimiGenerator.h"

namespace NetworKit {


 TEST_F(ConnectedComponentsGTest, testConnectedComponentsTiny) {
 	// construct graph
 	Graph g;
 	for (count i = 0; i < 20; i++) {
 		g.addNode();
 	}
 	g.addEdge(0,1,0);
 	g.addEdge(1,2,0);
 	g.addEdge(2,4,0);
 	g.addEdge(4,8,0);
 	g.addEdge(8,16,0);
 	g.addEdge(16,19,0);

 	g.addEdge(3,5,0);
 	g.addEdge(5,6,0);
 	g.addEdge(6,7,0);
 	g.addEdge(7,9,0);

 	g.addEdge(10,11,0);
 	g.addEdge(10,18,0);
 	g.addEdge(10,12,0);
 	g.addEdge(18,17,0);

 	g.addEdge(13,14,0);

 	// initialize ConnectedComponents
 	ConnectedComponents ccs(g);
 	ccs.run();

 	// check result
 	EXPECT_EQ(5, ccs.numberOfComponents());
 	EXPECT_TRUE(ccs.componentOfNode(0) == ccs.componentOfNode(19));
 	EXPECT_TRUE(ccs.componentOfNode(3) == ccs.componentOfNode(7));
 }


// TEST_F(ConnectedComponentsGTest, testGraphWithoutEdges) {
// 	// construct graph
// 	Graph g;
// 	for (count i = 0; i < 20; i++) {
// 		g.addNode();
// 	}

// 	// initialize ConnectedComponents
// 	ConnectedComponents ccs;
// 	ccs.run(g);

// 	// check result
// 	EXPECT_TRUE(ccs.numberOfComponents() == 20);
// 	EXPECT_TRUE(ccs.componentOfNode(0) != ccs.componentOfNode(19));
// 	EXPECT_TRUE(ccs.componentOfNode(3) != ccs.componentOfNode(7));
// }

// TEST_F(ConnectedComponentsGTest, testGetComponent) {
// 	// construct graph
// 	Graph g;
// 	for (count i = 0; i < 20; i++) {
// 		g.addNode();
// 	}
// 	g.addEdge(0,1,0);
// 	g.addEdge(1,2,0);
// 	g.addEdge(2,4,0);
// 	g.addEdge(4,8,0);
// 	g.addEdge(8,16,0);
// 	g.addEdge(16,19,0);

// 	g.addEdge(3,5,0);
// 	g.addEdge(5,6,0);
// 	g.addEdge(6,7,0);
// 	g.addEdge(7,9,0);

// 	g.addEdge(10,11,0);
// 	g.addEdge(10,18,0);
// 	g.addEdge(10,12,0);
// 	g.addEdge(18,17,0);

// 	g.addEdge(13,14,0);

// 	// initialize ConnectedComponents
// 	ConnectedComponents ccs;
// 	ccs.run(g);

// 	// check result
// 	std::vector<node> comp = ccs.getComponent(0);
// 	EXPECT_TRUE(comp.size() == ccs.sizeOfComponent(0));
// 	EXPECT_TRUE(comp.at(0) == 0);
// 	EXPECT_TRUE(comp.at(1) == 1);
// 	EXPECT_TRUE(comp.at(2) == 2);
// 	EXPECT_TRUE(comp.at(3) == 4);
// 	EXPECT_TRUE(comp.at(4) == 8);
// 	EXPECT_TRUE(comp.at(5) == 16);
// 	EXPECT_TRUE(comp.at(6) == 19);
// }

TEST_F(ConnectedComponentsGTest, testConnectedComponents) {
	// construct graph
	METISGraphReader reader;
	Graph G = reader.read("input/astro-ph.graph");
	ConnectedComponents cc(G);
	cc.run();
	DEBUG("Number of components: ", cc.numberOfComponents());
	EXPECT_EQ(1029u, cc.numberOfComponents());
}

TEST_F(ConnectedComponentsGTest, tryParallelConnectedComponents) {
	METISGraphReader reader;
	std::vector<std::string> graphs = {"astro-ph", "PGPgiantcompo",
			"caidaRouterLevel", "celegans_metabolic", "hep-th", "jazz"};

	for (auto graphName: graphs) {
		Graph G = reader.read("input/" + graphName + ".graph");
		ParallelConnectedComponents cc(G);
		cc.runSequential();
		count seqNum = cc.numberOfComponents();
		cc.run();
		count parNum = cc.numberOfComponents();
		DEBUG("Number of components: ", seqNum);
		EXPECT_EQ(seqNum, parNum);
	}
}

TEST_F(ConnectedComponentsGTest, benchConnectedComponents) {
	// construct graph
	METISGraphReader reader;
	Graph G = reader.read("input/coAuthorsDBLP.graph");
	ConnectedComponents cc(G);
	cc.run();
	DEBUG("Number of components: ", cc.numberOfComponents());
	EXPECT_EQ(1u, cc.numberOfComponents());
}

TEST_F(ConnectedComponentsGTest, benchHHConnectedComponents) {
	// construct graph
	METISGraphReader reader;
	Graph G = reader.read("input/coAuthorsDBLP.graph");
	ConnectedComponents cc(G);

	// run CC algo on original
	cc.run();
	DEBUG("Number of components in original: ", cc.numberOfComponents());

	// compute HH graph and apply CC algo
	std::vector<unsigned int> sequence = GraphProperties::degreeSequence(G);
	HavelHakimiGenerator hhgen(sequence, true);
	Graph G2 = hhgen.generate();
	ConnectedComponents cc2(G2);
	cc2.run();
	DEBUG("Number of components in HH generated: ", cc2.numberOfComponents());
}

TEST_F(ConnectedComponentsGTest, benchLiveJConnectedComponents) {
	// construct graph
	METISGraphReader reader;
	Graph G = reader.read("../graphs/ascii/soc-LiveJournal1_symm_or.graph");
	ConnectedComponents cc(G);

	// run CC algo on original
	cc.run();
	DEBUG("Number of components in original: ", cc.numberOfComponents());

	// compute HH graph and apply CC algo
	std::vector<unsigned int> sequence = GraphProperties::degreeSequence(G);
	HavelHakimiGenerator hhgen(sequence, true);
	Graph G2 = hhgen.generate();
	ConnectedComponents cc2(G2);
	cc2.run();
	DEBUG("Number of components in HH generated: ", cc2.numberOfComponents());
}


TEST_F(ConnectedComponentsGTest, testStronglyConnectedComponents) {

    auto comparePartitions = [](const Partition& p1, const Partition& p2) {
    std::vector<index> partitionIdMap(p1.upperBound(), none);
    ASSERT_EQ(p1.numberOfElements(), p2.numberOfElements());
    ASSERT_EQ(p1.numberOfSubsets(), p2.numberOfSubsets());

    p1.forEntries([&](node v, index p) {
        if (partitionIdMap[p] == none) {
            partitionIdMap[p] = p2.subsetOf(v);
        }
        index p_mapped = partitionIdMap[p];
        ASSERT_EQ(p_mapped, p);
    });
    };


    count n = 8;
    count m = 14;
    Graph G(n, false, true);

    G.addEdge(0, 4);
    G.addEdge(1, 0);
    G.addEdge(2, 1);
    G.addEdge(2, 3);
    G.addEdge(3, 2);
    G.addEdge(4, 1);
    G.addEdge(5, 1);
    G.addEdge(5, 4);
    G.addEdge(5, 6);
    G.addEdge(6, 2);
    G.addEdge(6, 5);
    G.addEdge(7, 3);
    G.addEdge(7, 6);
    G.addEdge(7, 7);

    ASSERT_EQ(n, G.numberOfNodes());
    ASSERT_EQ(m, G.numberOfEdges());

    count z = G.upperNodeIdBound();
    Partition p_expected(z);
    p_expected.allToSingletons();
    p_expected[0] = 0;
    p_expected[1] = 0;
    p_expected[2] = 1;
    p_expected[3] = 1;
    p_expected[4] = 0;
    p_expected[5] = 2;
    p_expected[6] = 2;
    p_expected[7] = 3;
    p_expected.compact();

    StronglyConnectedComponents scc(G);
    scc.run();
    Partition p_actual = scc.getPartition();
    p_actual.compact();

    comparePartitions(p_expected, p_actual);
}


} /* namespace NetworKit */

#endif /*NOGTEST */
