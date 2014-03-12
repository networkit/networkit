/*
 * ConnectedComponentsGTest.cpp
 *
 *  Created on: Sep 16, 2013
 *      Author: Maximilian Vogel
 */
#ifndef NOGTEST

#include "ConnectedComponentsGTest.h"
#include "../ConnectedComponents.h"
#include "../GraphProperties.h"
#include "../Diameter.h"
#include "../../io/METISGraphReader.h"
#include "../../generators/HavelHakimiGenerator.h"

namespace NetworKit {

ConnectedComponentsGTest::ConnectedComponentsGTest() {

}

ConnectedComponentsGTest::~ConnectedComponentsGTest() {

}

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
 	ConnectedComponents ccs;
 	ccs.run(g);

 	// check result
 	EXPECT_TRUE(ccs.numberOfComponents() == 5);
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

TEST_F(ConnectedComponentsGTest, testParallelConnectedComponents) {
	// construct graph
	METISGraphReader reader;
	Graph G = reader.read("input/PGPgiantcompo.graph");
	ConnectedComponents cc;
	cc.run(G);
	DEBUG("Number of components: ", cc.numberOfComponents());
	EXPECT_EQ(1, cc.numberOfComponents());
}

TEST_F(ConnectedComponentsGTest, tryHHConnectedComponents) {
	// construct graph
	METISGraphReader reader;
	Graph G = reader.read("input/coAuthorsDBLP.graph");
	ConnectedComponents cc;

	// run CC algo on original
	cc.run(G);
	DEBUG("Number of components in original: ", cc.numberOfComponents());

	// compute HH graph and apply CC algo
	std::vector<unsigned int> sequence = GraphProperties::degreeSequence(G);
	HavelHakimiGenerator hhgen(sequence, true);
	Graph G2 = hhgen.generate();
	cc.run(G2);
	DEBUG("Number of components in HH generated: ", cc.numberOfComponents());
}

TEST_F(ConnectedComponentsGTest, tryLiveJConnectedComponents) {
	// construct graph
	METISGraphReader reader;
	Graph G = reader.read("../graphs/ascii/soc-LiveJournal1_symm_or.graph");
	ConnectedComponents cc;

	// run CC algo on original
	cc.run(G);
	DEBUG("Number of components in original: ", cc.numberOfComponents());

	// compute HH graph and apply CC algo
	std::vector<unsigned int> sequence = GraphProperties::degreeSequence(G);
	HavelHakimiGenerator hhgen(sequence, true);
	Graph G2 = hhgen.generate();
	cc.run(G2);
	DEBUG("Number of components in HH generated: ", cc.numberOfComponents());
}


} /* namespace NetworKit */

#endif /*NOGTEST */

