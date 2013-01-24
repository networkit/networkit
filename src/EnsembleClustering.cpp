//============================================================================
// Name        : EnsembleClustering.cpp
// Author      : Christian Staudt
// Version     :
// Copyright   : © 2012, Christian Staudt
// Description : Hello World in C++, Ansi-style
//============================================================================

// includes
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
#include "aux/Timer.h"
#include "aux/Functions.h"
#include "graph/Graph.h"
#include "graph/GraphGenerator.h"
#include "ensemble/EnsembleClusterer.h"
#include "clustering/algo/LabelPropagation.h"
#include "clustering/base/Clustering.h"
#include "clustering/base/Modularity.h"
#include "clustering/base/ClusteringGenerator.h"
#include "io/METISGraphReader.h"





using namespace EnsembleClustering;




bool startWithGraph(Graph& G, int ensembleSize) {

	EnsembleClusterer ensemble;

		// CONFIGURE ENSEMBLE CLUSTERER

		// 1. Quality Measure
		QualityMeasure* qm = new Modularity();
		ensemble.setQualityMeasure(*qm);

		// 2. Base Clusterers
		for (int i = 0; i < ensembleSize; i += 1) {
			Clusterer* base = new LabelPropagation();
			ensemble.addBaseClusterer(*base);
		}

		// 3. Final Clusterer
		Clusterer* final = new LabelPropagation();
		ensemble.setFinalClusterer(*final);

		// RUN ENSEMBLE CLUSTERER
		Aux::Timer runtime;
		runtime.start();

		Clustering result = ensemble.run(G);

		runtime.stop();
		std::cout << "EnsembleClusterer runtime: " << runtime.elapsed().count() << " ms" << std::endl;


		// ANALYZE RESULT

		return true;
}

bool startWithPath(std::string graphPath, int ensembleSize) {
	assert (ensembleSize > 0);
	assert (! graphPath.empty());

	// READ GRAPH

	METISGraphReader reader;	// TODO: add support for multiple graph file formats

	// TIMING
	Aux::Timer readTimer;
	readTimer.start();
	//
	std::cout << "opening file... : " << graphPath << std::endl;

	Graph G = reader.read(graphPath);
	//
	readTimer.stop();
	std::cout << "read graph file in " << readTimer.elapsed().count() << " ms " << std::endl;
	// TIMING

	// startWithGraph(G, ensembleSize);
	return true;
}


bool startWithGenerated(int64_t n, int64_t k, double pin, double pout, int ensembleSize) {

	// prepare generated clustered random graph (planted partition)
	std::cout << "[BEGIN] making random clustering..." << std::flush;
	Graph emptyG(n);	// clustering generator needs a dummy graph
	ClusteringGenerator clusteringGen;
	Clustering planted = clusteringGen.makeRandomClustering(emptyG, k);
	std::cout << "[DONE]" << std::endl;

	std::cout << "[BEGIN] making clustered random graph..." << std::flush;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(planted, pin, pout);
	std::cout << "[DONE]" << std::endl;

	// call ensemble clusterer
	return startWithGraph(G, ensembleSize);
}



/**
 * Call this first to configure logging output.
 */
void configureLogging(std::string loglevel = "DEBUG") {
	// configure logging
	log4cxx::BasicConfigurator::configure();
	if (loglevel == "TRACE") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getTrace());

	} else if (loglevel == "DEBUG") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getDebug());
	} else if (loglevel == "INFO") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getInfo());
	} else if (loglevel == "WARN") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getWarn());
	} else if (loglevel == "ERROR") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getError());
	} else {
		std::cout << "unknown log level: " << loglevel;
		std::exit(1);
	}
}


// *** Option Parser Configuration ***//

class Arg: public OptionParser::Arg {

static OptionParser::ArgStatus Required(const OptionParser::Option& option, bool msg)
{
  if (option.arg != 0)
    return OptionParser::ARG_OK;

  if (msg) {
	  std::cout << "Option '" << option << "' requires an argument" << std::endl;
  }
  return OptionParser::ARG_ILLEGAL;
}

};


enum  optionIndex { UNKNOWN, HELP, LOGLEVEL, TESTS, GRAPH, GENERATE, ENSEMBLE_SIZE};
const OptionParser::Descriptor usage[] =
{
 {UNKNOWN, 0,"" , ""    ,OptionParser::Arg::None, "USAGE: EnsembleClustering [options]\n\n"
                                            "Options:" },
 {HELP,    0,"h" , "help",OptionParser::Arg::None, "  --help  \t Print usage and exit." },
 {LOGLEVEL,    0, "" , "loglevel", OptionParser::Arg::Required, "  --loglevel  \t set the log level" },
 {TESTS, 0, "t", "tests", OptionParser::Arg::None, "  --tests \t Run unit tests"},
 {GRAPH, 0, "g", "graph", OptionParser::Arg::Required, "  --graph \t Run ensemble clusterer on graph"},
 {GENERATE, 0, "", "generate", OptionParser::Arg::Required, "  --generate \t Run ensemble clusterer on generated graph with planted partition"},
 {ENSEMBLE_SIZE, 0, "", "ensemble-size", OptionParser::Arg::Required, "  --ensemble-size \t number of clusterers in the ensemble"},
 {UNKNOWN, 0,"" ,  ""   ,OptionParser::Arg::None, "\nExamples:\n"
                                            " TODO" },
 {0,0,0,0,0,0}
};





int main(int argc, char **argv) {
	std::cout << "*** EnsembleClustering: combining parallel clustering algorithms with an ensemble learning strategy *** " << std::endl;

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

	if (options[LOGLEVEL]) {
		configureLogging(options[LOGLEVEL].arg);
	} else {
		configureLogging();	// with default level
	}


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
		graphPath = options[GRAPH].arg;
		std::cout << "\t --graph=" << graphPath << std::endl;

		ensembleSize = std::atoi(options[ENSEMBLE_SIZE].arg);
		std::cout << "\t --ensemble-size=" << ensembleSize << std::endl;

		std::cout << "[START] EnsembleClustering --graph" << std::endl;

		bool done = startWithPath(graphPath, ensembleSize);

		if (done) {
		   std::cout << "[FINISH] EnsembleClustering --graph" << std::endl;
		}

	}

	if (options[GENERATE]) {
		// read n, k, pin, pout
		std::string genOption = options[GENERATE].arg;
		std::cout << "\t --generate=" << options[GENERATE].arg << std::endl;

		std::vector<std::string> genArgs = splitString(genOption, ',');
		assert (genArgs.size() == 4);
		int n = std::atoi(genArgs.at(0).c_str());
		int k = std::atoi(genArgs.at(1).c_str());
		double pin = std::atof(genArgs.at(2).c_str());
		double pout = std::atof(genArgs.at(3).c_str());


		ensembleSize = std::atoi(options[ENSEMBLE_SIZE].arg);
		std::cout << "\t --ensemble-size=" << ensembleSize << std::endl;

		bool done = startWithGenerated(n, k, pin, pout, ensembleSize);

		if (done) {
		   std::cout << "[FINISH] EnsembleClustering --generate" << std::endl;
		}
	}


}
