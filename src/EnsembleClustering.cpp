//============================================================================
// Name        : EnsembleClustering.cpp
// Author      : Christian Staudt
// Version     :
// Copyright   : © 2012, Christian Staudt
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <utility>
#include <unordered_map>

// log4cxx
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"



// GoogleTest
#include <gtest/gtest.h>

// EnsembleClustering
#include "aux/Log.h"
#include "aux/Noise.h"
#include "graph/Graph.h"
#include "io/METISParser.h"
#include "io/METISToGraph.h"
#include "matching/Matching.h"
#include "clustering/base/Clustering.h"
#include "clustering/base/ClusteringGenerator.h"
#include "graph/GraphGenerator.h"
#include "clustering/base/Modularity.h"

extern "C" {
#include "stinger.h"
}




using namespace EnsembleClustering;




/**
 * Call this first to configure logging output.
 */
void configureLogging() {
	// configure logging
	log4cxx::BasicConfigurator::configure();
	log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getDebug());
}




int main(int argc, char **argv) {

	std::cout << "running EnsembleClustering" << std::endl;

	configureLogging();

	INFO("=== starting unit tests ===");

	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();

}
