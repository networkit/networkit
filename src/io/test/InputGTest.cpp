/*
 * InputGTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: cls
 */

#include "InputGTest.h"

namespace EnsembleClustering {

TEST_F(InputGTest, testGraphIOEdgeList) {
	GraphGenerator graphGen;
	Graph G = graphGen.makeCircularGraph(20);
	GraphIO graphio;
	std::string path = "sandbox/edgelist.txt";
	graphio.toEdgeList(G, path);

	bool exists = false;
	std::ifstream file(path);
	if (file) {
		exists = true;
	}
	EXPECT_TRUE(exists) << "A file should have been created : " << path;
}


} /* namespace EnsembleClustering */
