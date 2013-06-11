/*
 * IOBenchmark.cpp
 *
 *  Created on: 01.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "IOBenchmark.h"

namespace NetworKit {

IOBenchmark::IOBenchmark() {
	// TODO Auto-generated constructor stub

}

IOBenchmark::~IOBenchmark() {
	// TODO Auto-generated destructor stub
}


TEST_F(IOBenchmark, timeMETISGraphReader) {
	std::string path = "";

	std::cout << "[INPUT] .graph file path >" << std::endl;
	std::getline(std::cin, path);

	Aux::Timer runtime;

	INFO("[BEGIN] reading graph: " << path);
	runtime.start();
	METISGraphReader reader;
	Graph G = reader.read(path);
	runtime.stop();
	INFO("[DONE] reading graph " << runtime.elapsedTag());

	EXPECT_TRUE(! G.isEmpty());

}

TEST_F(IOBenchmark, tryDGSReader) {
	// read example graph
	DGSReader reader;
	Graph G;
	GraphEventProxy Gproxy(G);
	reader.read("input/AuthorsGraph.dgs", Gproxy);

	// get input parameters
	count nodeCount = G.numberOfNodes();
	DEBUG("Number of nodes " << nodeCount);

}

} /* namespace NetworKit */
