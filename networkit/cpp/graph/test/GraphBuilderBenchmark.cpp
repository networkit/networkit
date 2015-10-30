/*
 * GraphBuilderBenchmark.cpp
 *
 *  Created on: 04.12.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#include "GraphBuilderBenchmark.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {


GraphBuilderBenchmark::GraphBuilderBenchmark() {
}

TEST_F(GraphBuilderBenchmark, benchmarkMETISReader) {
	METISGraphReader reader;
	measureInMs([&]() {
		auto G = reader.read("../algoDaten/graphs/eu-2005.graph");
		return G.numberOfNodes();
	}, 20);
}

} /* namespace NetworKit */

#endif /*NOGTEST */
