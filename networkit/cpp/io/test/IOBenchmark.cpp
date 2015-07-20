/*
 * IOBenchmark.cpp
 *
 *  Created on: 01.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "IOBenchmark.h"
#include "../RasterReader.h"
#include "../../generators/Quadtree/Quadtree.h"

namespace NetworKit {



TEST_F(IOBenchmark, timeMETISGraphReader) {
	std::string path = "";

	std::cout << "[INPUT] .graph file path >" << std::endl;
	std::getline(std::cin, path);

	Aux::Timer runtime;

	INFO("[BEGIN] reading graph: " , path);
	runtime.start();
	METISGraphReader reader;
	Graph G = reader.read(path);
	runtime.stop();
	INFO("[DONE] reading graph " , runtime.elapsedTag());

	EXPECT_TRUE(! G.isEmpty());

}


TEST_F(IOBenchmark, benchRasterReader) {
	double normalizationFactor = 0.05;
	RasterReader reader(normalizationFactor);
	std::vector<double> xcoords;
	std::vector<double> ycoords;
	Aux::Timer runtime;
	std::string path("input/deup00ag.asc");

	// read raster file
	INFO("[BEGIN] reading raster data set: ", path);
	runtime.start();
	std::tie(xcoords, ycoords) = reader.read(path);
	runtime.stop();
	INFO("[DONE] reading raster data set " , runtime.elapsedTag());
	EXPECT_EQ(xcoords.size(), ycoords.size());

	// perform range queries
	// Quadtree(xcoords, ycoords);
	uint64_t numQueries = 1000;
	for (uint64_t q = 0; q < numQueries; ++q) {
		// TODO
	}
}


} /* namespace NetworKit */


#endif /* NOGTEST */
