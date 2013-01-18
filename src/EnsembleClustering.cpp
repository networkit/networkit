//============================================================================
// Name        : EnsembleClustering.cpp
// Author      : Christian Staudt
// Version     :
// Copyright   : © 2012, Christian Staudt
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <utility>


// log4cxx
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"

// GoogleTest
#include <gtest/gtest.h>


// EnsembleClustering
#include "aux/optionparser.h"
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




using namespace EnsembleClustering;


/**
 * Call this first to configure logging output.
 */
void configureLogging() {
	// configure logging
	log4cxx::BasicConfigurator::configure();
	log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getDebug());	// TODO: make debug level a command line parameter
}


// *** Option Parser Configuration ***//

class Arg: public OptionParser::Arg {

static OptionParser::ArgStatus Required(const OptionParser::Option& option, bool msg)
{
  if (option.arg != 0)
    return OptionParser::ARG_OK;

  if (msg) {
	  ERROR("Option '" << option << "' requires an argument\n");
  }
  return OptionParser::ARG_ILLEGAL;
}

};


enum  optionIndex { UNKNOWN, HELP, TESTS, GRAPH, ENSEMBLE_SIZE};
const OptionParser::Descriptor usage[] =
{
 {UNKNOWN, 0,"" , ""    ,OptionParser::Arg::None, "USAGE: EnsembleClustering [options]\n\n"
                                            "Options:" },
 {HELP,    0,"h" , "help",OptionParser::Arg::None, "  --help  \t Print usage and exit." },
 {TESTS, 0, "t", "tests", OptionParser::Arg::None, "  --tests \t Run unit tests"},
 {GRAPH, 0, "g", "graph", OptionParser::Arg::Required, "  --graph \t Run ensemble clusterer on graph"},
 {ENSEMBLE_SIZE, 0, "", "ensemble-size", OptionParser::Arg::Required, "  --ensemble-size \t number of clusterers in the ensemble"},
 {UNKNOWN, 0,"" ,  ""   ,OptionParser::Arg::None, "\nExamples:\n"
                                            " TODO" },
 {0,0,0,0,0,0}
};


void startEnsembleClusterer(std::string graphPath, int ensembleSize) {
	assert (ensembleSize > 0);
}



int main(int argc, char **argv) {
	/// PARSE OPTIONS

	argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
	OptionParser::Stats  stats(usage, argc, argv);
	OptionParser::Option options[stats.options_max], buffer[stats.buffer_max];
	OptionParser::Parser parse(usage, argc, argv, options, buffer);

	if (parse.error())
	 return 1;

	if (options[HELP] || argc == 0) {
	 OptionParser::printUsage(std::cout, usage);
	 return 0;
	}

	for (OptionParser::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
	 std::cout << "Unknown option: " << opt->name << "\n";

	for (int i = 0; i < parse.nonOptionsCount(); ++i)
	 std::cout << "Non-option #" << i << ": " << parse.nonOption(i) << "\n";


	// CONFIGURE LOGGING
	configureLogging();

	// RUN UNIT TESTS

	if (options[TESTS]) {

	   ::testing::InitGoogleTest(&argc, argv);
	   INFO("=== starting unit tests ===");
	   return RUN_ALL_TESTS();
	}

	// RUN PROGRAM ON GRAPH INSTANCE

	int ensembleSize = -1;
	std::string graphPath = "NONE";

	if (options[GRAPH]) {
	   DEBUG("called with graph argument: " << options[GRAPH].arg);
	}

	if (options[ENSEMBLE_SIZE]) {
	   ensembleSize = std::atoi(options[ENSEMBLE_SIZE].arg);
	   DEBUG("called with ensemble size: " << ensembleSize);
	}

	if ((graphPath != "") || (ensembleSize > 0)) {
	   INFO("starting EnsembleClusterer");

	} else {
	   ERROR("wrong options");
	   return 1;
	}


}
