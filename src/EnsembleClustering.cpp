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

// import global



using namespace EnsembleClustering;



// graph functions

Graph generateGraph(int64_t n, int64_t k, double pin, double pout) {

	// prepare generated clustered random graph (planted partition)
	std::cout << "[BEGIN] making random clustering..." << std::flush;
	Graph emptyG(n);	// clustering generator needs a dummy graph
	ClusteringGenerator clusteringGen;
	Clustering planted = clusteringGen.makeRandomClustering(emptyG, k);
	std::cout << "[DONE]" << std::endl;

	std::cout << "[BEGIN] making clustered random graph..." << std::flush;
	GraphGenerator graphGen;

	Aux::Timer runtime;
	runtime.start();
	Graph G = graphGen.makeClusteredRandomGraph(planted, pin, pout);
	runtime.stop();
	std::cout << "[DONE] (" << runtime.elapsed().count() << " ms)" << std::endl;

	return G;
}


Graph readGraph(std::string graphPath) {

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

	return G;

}

// start functions
// FIXME: deprecated
bool startSoloWithGraph(Graph& G, Clusterer* algo) {

	// RUN BASE CLUSTERER SOLO
	Aux::Timer runtime;
	runtime.start();

	Clustering zeta = algo->run(G);

	runtime.stop();
	std::cout << "Solo clusterer runtime: " << runtime.elapsed().count() << " ms" << std::endl;

	return true;
}


// FIXME: deprecated
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

// FIXME: deprecated
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

// FIXME: deprecated
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




std::vector<std::string> splitAString(std::string s, char delim = ' ') {
	std::stringstream stream(s);
	std::string token;
	std::vector<std::string> tokens;

	// split string and push adjacent nodes
	while (std::getline(stream, token, delim)) {
		tokens.push_back(token);
	}

	return tokens;
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


enum  optionIndex { UNKNOWN, HELP, LOGLEVEL, SEQUENTIAL, TESTS, GRAPH, GENERATE, ENSEMBLE_SIZE, ENSEMBLE, SOLO};
const OptionParser::Descriptor usage[] =
{
 {UNKNOWN, 0,"" , ""    ,OptionParser::Arg::None, "USAGE: EnsembleClustering [options]\n\n"
                                            "Options:" },
 {HELP,    0,"h" , "help",OptionParser::Arg::None, "  --help  \t Print usage and exit." },
 {LOGLEVEL,    0, "" , "loglevel", OptionParser::Arg::Required, "  --loglevel  \t set the log level" },
 {SEQUENTIAL,    0, "s" , "sequential", OptionParser::Arg::None, "  --sequential  \t Turn off parallelism." },
 {TESTS, 0, "t", "tests", OptionParser::Arg::None, "  --tests \t Run unit tests"},
 {GRAPH, 0, "g", "graph", OptionParser::Arg::Required, "  --graph \t Run ensemble clusterer on graph"},
 {GENERATE, 0, "", "generate", OptionParser::Arg::Required, "  --generate \t Run ensemble clusterer on generated graph with planted partition"},
 {ENSEMBLE_SIZE, 0, "", "ensemble-size", OptionParser::Arg::Required, "  --ensemble-size \t number of clusterers in the ensemble"},
 {ENSEMBLE, 0, "", "ensemble", OptionParser::Arg::Required, "  --ensemble=<b> \t <b>: number of base clusterers in the ensemble"},	// TODO: provide more options
 {SOLO, 0, "", "solo", OptionParser::Arg::Required, "  --solo=<Algorithm> \t run only a single base algorithm"},
 {UNKNOWN, 0,"" ,  ""   ,OptionParser::Arg::None, "\nExamples:\n"
                                            " TODO" },
 {0,0,0,0,0,0}
};


// MAIN FUNCTIONS


Graph getGraph(OptionParser::Option* options) {

	// generated graph
	if (options[GENERATE]) {
		// read n, k, pin, pout
		std::string genOption = options[GENERATE].arg;
		std::cout << "\t --generate=" << options[GENERATE].arg << std::endl;

		std::vector<std::string> genArgs = splitAString(genOption, ',');
		assert (genArgs.size() == 4);
		int64_t n = std::atoi(genArgs.at(0).c_str());
		int64_t k = std::atoi(genArgs.at(1).c_str());
		double pin = std::atof(genArgs.at(2).c_str());
		double pout = std::atof(genArgs.at(3).c_str());

		Graph G = generateGraph(n, k, pin, pout);
		return G;
	} else if (options[GRAPH]) { // graph from file
		std::string graphPath = options[GRAPH].arg;
		std::cout << "\t --graph=" << graphPath << std::endl;

		Graph G = readGraph(graphPath);
		return G;
	} else {
		Graph G(0);	// return empty graph
		G.setName("NONE");
		return G;
	}
}


bool startAlgo(Graph G, OptionParser::Option* options) {

	// if getGraph returns empty graph, abort
	if (G.isEmpty() && (G.getName() == "NONE")) {
		std::cout << "[ERROR]Êno graph instance provided." << std::endl;
		std::cout << "[EXIT]" << std::endl;
		exit(1);
	}


	if (options[SOLO]) {
		std::cout << "\t --solo=" << options[SOLO].arg << std::endl;
		// RUN ONLY SINGLE BASE ALGORITHM

		Clusterer* algo = NULL; // the clusterer

		// get specified base algorithm
		std::string algoName = options[SOLO].arg;
		if (algoName == "LabelPropagation") {
			algo = new LabelPropagation();
		} else {
			std::cout << "[ERROR] unknown base algorithm: " << algoName << std::endl;
			std::cout << "[EXIT]" << std::endl;
			exit(1);
		}

		// start solo base algorithm
		std::cout << "[BEGIN] solo base clusterer: " << algoName << std::endl;

		algo->run(G);

		//
		std::cout << "[DONE]" << std::endl;

	} else if (options[ENSEMBLE]) {
		// RUN ENSEMBLE

		std::string ensembleOptions = options[GENERATE].arg;

		int ensembleSize = std::atoi(ensembleOptions.c_str()); // TODO: provide more options

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


	}

}



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


	// PARALLELISMISM SWITCH
	if (options[SEQUENTIAL]) {
		// TODO: parallelism switch
	}


	// RUN UNIT TESTS

	if (options[TESTS]) {

	   ::testing::InitGoogleTest(&argc, argv);
	   INFO("=== starting unit tests ===");
	   return RUN_ALL_TESTS();
	}

	//




	// RUN PROGRAM

	startAlgo(getGraph(options), options);


//	int ensembleSize = -1;
//	std::string graphPath = "NONE";
//
//	if (options[GRAPH]) {
//		graphPath = options[GRAPH].arg;
//		std::cout << "\t --graph=" << graphPath << std::endl;
//
//		ensembleSize = std::atoi(options[ENSEMBLE_SIZE].arg);
//		std::cout << "\t --ensemble-size=" << ensembleSize << std::endl;
//
//		std::cout << "[START] EnsembleClustering --graph" << std::endl;
//
//		bool done = startWithPath(graphPath, ensembleSize);
//
//		if (done) {
//		   std::cout << "[FINISH] EnsembleClustering --graph" << std::endl;
//		}
//	}
//
//	if (options[GENERATE]) {
//		// read n, k, pin, pout
//		std::string genOption = options[GENERATE].arg;
//		std::cout << "\t --generate=" << options[GENERATE].arg << std::endl;
//
//		std::vector<std::string> genArgs = splitAString(genOption, ',');
//		assert (genArgs.size() == 4);
//		int n = std::atoi(genArgs.at(0).c_str());
//		int k = std::atoi(genArgs.at(1).c_str());
//		double pin = std::atof(genArgs.at(2).c_str());
//		double pout = std::atof(genArgs.at(3).c_str());
//
//
//		ensembleSize = std::atoi(options[ENSEMBLE_SIZE].arg);
//		std::cout << "\t --ensemble-size=" << ensembleSize << std::endl;
//
//		bool done = startWithGenerated(n, k, pin, pout, ensembleSize);
//
//		if (done) {
//		   std::cout << "[FINISH] EnsembleClustering --generate" << std::endl;
//		}
//	}
//
//
//	// RUN ONLY SINGLE BASE ALGORITHM
//	if (options[SOLO]) {
//		Clusterer* algo = NULL;
//
//		std::string algoName = options[SOLO].arg;
//		if (algoName == "LabelPropagation") {
//			algo = new LabelPropagation();
//		} else {
//			std::cout << "[ERROR] unknown base algorithm: " << algoName << std::endl;
//			std::cout << "[EXIT]" << std::endl;
//			exit(1);
//		}
//
//		if (options[GRAPH]) {
//			graphPath = options[GRAPH].arg;
//			// TODO: start solo with graph
//		}
//
//	}


}
