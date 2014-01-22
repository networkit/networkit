/*
 * ConductanceTreeGTest.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Henning
 */

#ifndef NOGTEST

#include "ConductanceTreeGTest.h"


namespace NetworKit {


ConductanceTreeGTest::ConductanceTreeGTest() {

}

ConductanceTreeGTest::~ConductanceTreeGTest() {

}


TEST_F(ConductanceTreeGTest, tryConductanceTree) {
	// read application graph
	METISGraphReader graphReader;
	Graph graph = graphReader.read("input/PGPgiantcompo.graph");

	// check and evaluate mapping
	ConductanceTree tree;
	node root = 0;
	Clustering condBest = tree.bestCutInBfsTree(graph, root);
//	double condVal = ;
//	INFO("Conductance value of best cut found: " , condVal);

	// TODO evaluate conductance with method in clustering; test fixed value
}


} // namespace EnsembleClustering

#endif
