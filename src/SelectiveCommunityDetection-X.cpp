//============================================================================
// Name        : SelectiveCommunityDetection.cpp
// Author      : Christian Staudt (christian.staudt@kit.edu),
//				 Henning Meyerhenke (henning.meyerhenke@kit.edu)
// Version     :
// Copyright   : � 2013, Christian Staudt, Henning Meyerhenke
// Description : Framework for engineering selective community detection algorithms
//============================================================================

// includes

// includes
#include <iostream>
#include <fstream>
#include <utility>
#include <cfenv>	// floating point exceptions
#include <stdexcept>

// log4cxx
#ifndef NOLOG4CXX
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#endif

// GoogleTest
#ifndef NOGTEST
#include <gtest/gtest.h>
#endif

// OpenMP
#include <omp.h>

// EnsembleClustering
#include "Globals.h"
#include "ext/optionparser.h"
#include "auxiliary/Log.h"
#include "auxiliary/Timer.h"
#include "auxiliary/Functions.h"
#include "auxiliary/StringTools.h"
#include "graph/Graph.h"
#include "graph/Subgraph.h"
#include "io/METISGraphReader.h"
#include "io/METISGraphWriter.h"
#include "scd/RandomSeedSet.h"
#include "scd/RandomWalkSeedSet.h"
#include "scd/TSelectiveSCAN.h"
#include "scd/SelectiveSCAN.h"
#include "scd/TGreedyCommunityExpansion.h"
#include "scd/TQualityObjective.h"
#include "scd/TAcceptability.h"
#include "scd/GreedyCommunityExpansion.h"
#include "scd/QualityObjective.h"
#include "scd/Acceptability.h"
#include "scd/CommunityTrimming.h"
#include "distmeasures/TAlgebraicDistance.h"
#include "distmeasures/TNeighborhoodDistance.h"
#include "distmeasures/TNodeDistance.h"
#include "distmeasures/AlgebraicDistance.h"
#include "distmeasures/NeighborhoodDistance.h"
#include "distmeasures/NodeDistance.h"
#include "scd/CommunityQualityMeasure.h"

using namespace NetworKit;

#ifndef NOLOGGING
#ifndef NOLOG4CXX
/**
 * Call this first to configure logging output.
 */
void configureLogging(const std::string& loglevel = "INFO") {
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
		ERROR("unknown loglevel: " << loglevel);
		exit(1);
	}
}
#endif
#endif

/**
 *  Set the number of threads available to the program.
 */
void setNumberOfThreads(int nThreads) {
#ifdef _OPENMP
	omp_set_num_threads(nThreads);
#else
	WARN("Thread option ignored since OpenMP is deactivated.");
#endif
}

// *** Option Parser Configuration ***//

class Arg: public OptionParser::Arg {

	static OptionParser::ArgStatus Required(const OptionParser::Option& option,
			bool msg) {
		if (option.arg != 0)
			return OptionParser::ARG_OK;

		if (msg) {
			std::cout << "Option '" << option << "' requires an argument"
					<< std::endl;
		}
		return OptionParser::ARG_ILLEGAL;
	}

};

// TODO: clean up obsolete parameters
enum optionIndex {
	UNKNOWN,
	HELP,
	LOGLEVEL,
	THREADS,
	TESTS,
	GRAPH,
	SEEDS,
	DETECTOR,
	PARAM,
	QUALITY,
	GROUND_TRUTH,
	RUNS,
	SAVE_GRAPH,
	PROGRESS,
	SUMMARY,
	SHOW
};
const OptionParser::Descriptor usage[] =
		{ { UNKNOWN, 0, "", "", OptionParser::Arg::None,
				"USAGE: EnsembleClustering [options]\n\n"
						"Options:" }, { HELP, 0, "h", "help",
				OptionParser::Arg::None, "  --help  \t Print usage and exit." },
				{ LOGLEVEL, 0, "", "loglevel", OptionParser::Arg::Required,
						"  --loglevel=<LEVEL>  \t set the log level" },
				{ THREADS, 0, "", "threads", OptionParser::Arg::Required,
						"  --threads=<NUM>  \t set the maximum number of threads" },
				{ TESTS, 0, "t", "tests", OptionParser::Arg::None,
						"  --tests \t Run unit tests" }, { GRAPH, 0, "g",
						"graph", OptionParser::Arg::Required,
						"  --graph=<PATH> \t Run ensemble clusterer on graph" },
				{ SEEDS, 0, "", "seeds", OptionParser::Arg::Required,
						"  --seeds=<NAME> \t Specify seed set generator" },
				{ DETECTOR, 0, "", "detector", OptionParser::Arg::Required,
						"  --detector=<NAME>:<PARAMS> \t select clustering algorithm" },
				{ PARAM, 0, "", "param", OptionParser::Arg::Required,
						"  --param=<NAME>:<PARAMS> \t select parameters" },
				{ QUALITY, 0, "", "quality", OptionParser::Arg::Required,
						"  --quality=<NAME> \t select quality measure for evaluation" },
				{ GROUND_TRUTH, 0, "", "groundTruth",
						OptionParser::Arg::Required,
						"  --groundTruth=<NAME>:<PARAMS> \t select ground truth" },
				{ RUNS, 0, "", "runs", OptionParser::Arg::Required,
						"  --runs=<NUMBER> \t set number of clusterer runs" }, {
						SAVE_GRAPH, 0, "", "saveGraph",
						OptionParser::Arg::Required,
						"  --saveGraph=<PATH> \t write the graph to a file" }, {
						PROGRESS, 0, "", "progress", OptionParser::Arg::None,
						"  --progress \t print progress bar" },
				{ SUMMARY, 0, "", "summary", OptionParser::Arg::Required,
						"  --summary=<PATH> \t append summary as a .csv line to this file" },
				{ SHOW, 0, "", "show", OptionParser::Arg::Required,
						"  --show=<NAME>:<PARAMS> \t select parameters" },
				{ UNKNOWN, 0, "", "", OptionParser::Arg::None, "\nExamples:\n"
						" TODO" }, { 0, 0, 0, 0, 0, 0 } };

// MAIN FUNCTIONS

/**
 * Read a graph from a file.
 */
Graph readGraph(const std::string& graphPath) {

	// READ GRAPH

	METISGraphReader reader; // TODO: add support for multiple graph file formats

	// TIMING
	Aux::Timer readTimer;
	readTimer.start();
	//
	std::cout << "[BEGIN] reading file: " << graphPath << std::endl;

	Graph G = reader.read(graphPath);
	//
	readTimer.stop();
	std::cout << "[DONE] read graph file " << readTimer.elapsedTag()
			<< std::endl;
	// TIMING

	return G;

}

Graph getGraph(OptionParser::Option* options) {

	if (options[GRAPH]) { // graph from file
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

int main(int argc, char **argv) {
	std::cout << "=== NetworKit - Selective Community Detection ==="
			<< std::endl;

	// ENABLE FLOATING POINT EXCEPTIONS (needs GNU extension, apparently only available on Linux)
#ifdef _GNU_SOURCE
	// feenableexcept(FE_ALL_EXCEPT);
#endif

	/// PARSE OPTIONS

	argc -= (argc > 0);
	argv += (argc > 0); // skip program name argv[0] if present
	OptionParser::Stats stats(usage, argc, argv);
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

#ifndef NOLOGGING
#ifndef NOLOG4CXX
	if (options[LOGLEVEL]) {
		configureLogging(options[LOGLEVEL].arg);
	} else {
		configureLogging();	// with default level
	}
#endif
#endif

	// CONFIGURE PARALLELISM

#ifdef _OPENMP
	omp_set_nested(1); // enable nested parallelism
#endif

	if (options[THREADS]) {
		// set number of threads
		int nThreads = std::atoi(options[THREADS].arg);
		setNumberOfThreads(nThreads);
	}

	// CONFIGURE OUTPUT
	if (options[PROGRESS]) {
		PRINT_PROGRESS = true;
	} else {
		PRINT_PROGRESS = false;
	}


	// RUN UNIT TESTS

	if (options[TESTS]) {

#ifndef NOGTEST
		::testing::InitGoogleTest(&argc, argv);
		INFO("=== starting unit tests ===");
		return RUN_ALL_TESTS();
#else
		throw std::runtime_error("unit tests are excluded from build by the NOGTEST preprocessor directive");
#endif
	}

	// SET NUMBER OF CLUSTERER RUNS
	int runs = 1;
	if (options[RUNS]) {
		runs = std::atoi(options[RUNS].arg);
	}

	// RUN PROGRAM

	// get graph
	Graph G = getGraph(options);
	SeedSetGenerator* seedGen = NULL;
	count nSeeds = 1; // number of seeds

	// get seed set generator
	if (options[SEEDS]) {
		std::string seedsArg = options[SEEDS].arg;
		std::string seedsName = Aux::StringTools::split(seedsArg, ':')[0];
		std::string seedsParam = Aux::StringTools::split(seedsArg, ':')[1];

		nSeeds = std::stoi(seedsParam); // number of seeds

		if (seedsName == "RandomSeedSet") {
			seedGen = new RandomSeedSet(G);
		} else if (seedsName == "RandomWalkSeedSet") {
			seedGen = new RandomWalkSeedSet(G);
		} else {
			std::cout << "[ERROR] unknown argument for option --seeds: "
					<< seedsName << std::endl;
			exit(1);
		}

	} else {
		std::cout << "[ERROR] option --seeds must be supplied" << std::endl;
		exit(1);
	}



	Parameters param;
	double epsilon = 0.3;
	count mu = 2;
	double omega = 0.0;
	count numSystems = 0;
	count numIters = 0;
	count norm = 0;

	if (options[PARAM]) {
		std::string paramArg = options[PARAM].arg;
		std::vector<std::string> paramVec = Aux::StringTools::split(paramArg, ',');
		while (paramVec.size() > 0) {
			TRACE("parsing " << paramArg);
			std::string parameter = paramVec[0];
			if (Aux::StringTools::split(parameter, ':').size() == 2) {
				std::string first =
						Aux::StringTools::split(parameter, ':').front();
				std::string second =
						Aux::StringTools::split(parameter, ':').back();
				if (first == "epsilon") {
					epsilon = std::stof(second);
				} else if (first == "mu") {
					mu = std::stoi(second);
				} else if (first == "numSystems") {
					numSystems = std::stoi(second);
				} else if (first == "numIters") {
					numIters = std::stoi(second);
				} else if (first == "omega") {
					omega = std::stof(second);
				} else if (first == "norm") {
					norm = std::stoi(second);
				} else {
					std::cout << "[ERROR] invalid argument: " << second << std::endl;
					exit(1);
				}
			} else {
				std::cout << "[ERROR] invalid argument: " << parameter << std::endl;
				exit(1);
			}
			paramVec.erase(paramVec.begin());
		}
	}

	param.setInt("norm", norm);
	param.setInt("numSystems", numSystems);
	param.setInt("numIters", numIters);
	param.setDouble("omega", omega);


	std::vector<CommunityQualityMeasure*> measures;
	if (options[QUALITY]) {
		std::string qualityArg = options[QUALITY].arg;
		std::vector<std::string> qualityVec = Aux::StringTools::split(qualityArg, ',');
		while (qualityVec.size() > 0) {
			std::string parameter = qualityVec[0];
			if (parameter == "Conductance") {
				CommunityQualityMeasure* measure = NULL;
				measure = new Conduct(G);
				measures.push_back(measure);
			} else if ("Modularity") {
				CommunityQualityMeasure* measure = NULL;
				measure = new LocalModularity(G);
				measures.push_back(measure);
			} else {
				std::cout << "[ERROR] unknown quality measure: " << parameter << std::endl;
				exit(1);
			}
			qualityVec.erase(qualityVec.begin());
		}
	}

	// get algorithms
	SelectiveCommunityDetector* algo = NULL;
	std::unordered_set<node> tmp;
	std::unordered_map<node, count> bound;
	Acceptability* similarity = NULL;
	QualityObjective* objective = NULL;
	CommunityTrimming* trimming = NULL;
	NodeDistance* dist = NULL;

	if (options[DETECTOR]) {
		std::string detectorArg = options[DETECTOR].arg;
		std::string detectorName = Aux::StringTools::split(detectorArg, ':')[0];

		if (detectorName == "TSelectiveSCAN") {
			if (Aux::StringTools::split(detectorArg, ':').size() > 1) {
				if (Aux::StringTools::split(detectorArg, ':').back() == "TAD") {
					algo = new TSelectiveSCAN<TAlgebraicDistance>(G, param,
							epsilon, mu);
				} else {
					algo = new TSelectiveSCAN<TNeighborhoodDistance>(G, param,
							epsilon, mu);
				}
			} else {
				algo = new TSelectiveSCAN<TNeighborhoodDistance>(G, param,
						epsilon, mu);
			}
		} else if (detectorName == "SelectiveSCAN") {
			if (Aux::StringTools::split(detectorArg, ':').size() > 1) {
				if (Aux::StringTools::split(detectorArg, ':').back() == "AD") {
					dist = new AlgebraicDistance(G, numSystems, numIters, omega, norm);
					algo = new SelectiveSCAN(G, *dist, epsilon, mu);
				} else {
					dist = new NeighborhoodDistance(G);
					algo = new SelectiveSCAN(G, *dist, epsilon, mu);
				}
			} else {
				dist = new NeighborhoodDistance(G);
				algo = new SelectiveSCAN(G, *dist, epsilon, mu);
			}
		} else if (detectorName == "TGreedyCommunityExpansion") {
			if (Aux::StringTools::split(detectorArg, ':').size() == 4) {
				std::string first = Aux::StringTools::split(detectorArg, ':')[1];
				std::string second = Aux::StringTools::split(detectorArg, ':')[2];
				std::string third = Aux::StringTools::split(detectorArg, ':')[3];
				if (first == "Conductance" && second == "NodeClusterSimilarity"
						&& third == "BoundarySharpness") {
					algo = new TGreedyCommunityExpansion<TConductance,
							TNodeClusterSimilarity, BoundarySharpness>(G);
				} else if (first == "Conductance"
						&& second == "NodeClusterSimilarity"
						&& third == "Dummy") {
					algo = new TGreedyCommunityExpansion<TConductance,
							TNodeClusterSimilarity, DummyTrimming>(G);
				} else if (first == "Conductance" && second == "Dummy"
						&& third == "BoundarySharpness") {
					algo = new TGreedyCommunityExpansion<TConductance,
							TDummyAcceptability, BoundarySharpness>(G);
				} else if (first == "Conductance" && second == "Dummy"
						&& third == "Dummy") {
					algo = new TGreedyCommunityExpansion<TConductance,
							TDummyAcceptability, DummyTrimming>(G);
				} else if (first == "ModularityM"
						&& second == "NodeClusterSimilarity"
						&& third == "BoundarySharpness") {
					algo = new TGreedyCommunityExpansion<TLocalModularityM,
							TNodeClusterSimilarity, BoundarySharpness>(G);
				} else if (first == "ModularityM"
						&& second == "NodeClusterSimilarity"
						&& third == "Dummy") {
					algo = new TGreedyCommunityExpansion<TLocalModularityM,
							TNodeClusterSimilarity, DummyTrimming>(G);
				} else if (first == "ModularityM" && second == "Dummy"
						&& third == "BoundarySharpness") {
					algo = new TGreedyCommunityExpansion<TLocalModularityM,
							TDummyAcceptability, BoundarySharpness>(G);
				} else if (first == "ModularityM" && second == "Dummy"
						&& third == "Dummy") {
					algo = new TGreedyCommunityExpansion<TLocalModularityM,
							TDummyAcceptability, DummyTrimming>(G);
				} else if (first == "ModularityL"
						&& second == "NodeClusterSimilarity"
						&& third == "BoundarySharpness") {
					algo = new TGreedyCommunityExpansion<TLocalModularityL,
							TNodeClusterSimilarity, BoundarySharpness>(G);
				} else if (first == "ModularityL"
						&& second == "NodeClusterSimilarity"
						&& third == "Dummy") {
					algo = new TGreedyCommunityExpansion<TLocalModularityL,
							TNodeClusterSimilarity, DummyTrimming>(G);
				} else if (first == "ModularityL" && second == "Dummy"
						&& third == "BoundarySharpness") {
					algo = new TGreedyCommunityExpansion<TLocalModularityL,
							TDummyAcceptability, BoundarySharpness>(G);
				} else if (first == "ModularityL" && second == "Dummy"
						&& third == "Dummy") {
					algo = new TGreedyCommunityExpansion<TLocalModularityL,
							TDummyAcceptability, DummyTrimming>(G);
				} else {
					std::cout << "[ERROR] invalid arguments: " << first << " or " << second << " or " << third << std::endl;
					exit(1);
				}

			} else if (Aux::StringTools::split(detectorArg, ':').size() == 3) {
				std::string first = Aux::StringTools::split(detectorArg, ':')[1];
				std::string second =
						Aux::StringTools::split(detectorArg, ':')[2];
				if (first == "Conductance") {
					if (second == "NodeClusterSimilarity") {
						algo = new TGreedyCommunityExpansion<TConductance,
								TNodeClusterSimilarity, DummyTrimming>(G);
					} else if (second == "BoundarySharpness") {
						algo = new TGreedyCommunityExpansion<TConductance,
								TDummyAcceptability, BoundarySharpness>(G);
					} else if (second == "Dummy") {
						algo = new TGreedyCommunityExpansion<TConductance,
								TDummyAcceptability, DummyTrimming>(G);
					} else {
						std::cout << "[ERROR] invalid arguments " << std::endl;
						exit(1);
					}
				} else if (first == "ModularityM") {
					if (second == "NodeClusterSimilarity") {
						algo = new TGreedyCommunityExpansion<TLocalModularityM,
								TNodeClusterSimilarity, DummyTrimming>(G);
					} else if (second == "BoundarySharpness") {
						algo = new TGreedyCommunityExpansion<TLocalModularityM,
								TDummyAcceptability, BoundarySharpness>(G);
					} else if (second == "Dummy") {
						algo = new TGreedyCommunityExpansion<TLocalModularityM,
								TDummyAcceptability, DummyTrimming>(G);
					} else {
						std::cout << "[ERROR] invalid arguments " << std::endl;
						exit(1);
					}
				} else if (first == "ModularityL") {
					if (second == "NodeClusterSimilarity") {
						algo = new TGreedyCommunityExpansion<TLocalModularityL,
								TNodeClusterSimilarity, DummyTrimming>(G);
					} else if (second == "BoundarySharpness") {
						algo = new TGreedyCommunityExpansion<TLocalModularityL,
								TDummyAcceptability, BoundarySharpness>(G);
					} else if (second == "Dummy") {
						algo = new TGreedyCommunityExpansion<TLocalModularityL,
								TDummyAcceptability, DummyTrimming>(G);
					} else {
						std::cout << "[ERROR] invalid arguments " << std::endl;
						exit(1);
					}
				} else {
					std::cout << "[ERROR] invalid arguments: " << first << " or " << second << std::endl;
					exit(1);
				}
			} else if (Aux::StringTools::split(detectorArg, ':').size() == 2) {
				if (Aux::StringTools::split(detectorArg, ':')[1]
						== "Conductance") {
					algo = new TGreedyCommunityExpansion<TConductance,
							TDummyAcceptability, DummyTrimming>(G);
				} else if (Aux::StringTools::split(detectorArg, ':')[1]
						== "ModularityM") {
					algo = new TGreedyCommunityExpansion<TLocalModularityM,
							TDummyAcceptability, DummyTrimming>(G);
				} else if (Aux::StringTools::split(detectorArg, ':')[1]
						== "ModularityL") {
					algo = new TGreedyCommunityExpansion<TLocalModularityL,
							TDummyAcceptability, DummyTrimming>(G);
				} else {
					std::cout << "[ERROR] invalid arguments: " << detectorArg << std::endl;
					exit(1);
				}
			} else {
				std::cout << "[ERROR] invalid argument " << std::endl;
				exit(1);
			}
		} else if (detectorName == "GreedyCommunityExpansion") {

			if (Aux::StringTools::split(detectorArg, ':').size() == 4) {
				std::string first = Aux::StringTools::split(detectorArg, ':')[1];
				std::string second =
						Aux::StringTools::split(detectorArg, ':')[2];
				std::string third = Aux::StringTools::split(detectorArg, ':')[3];
				if (first == "Conductance" && second == "NodeClusterSimilarity"
						&& third == "BoundarySharpness") {
					similarity = new NodeClusterSimilarity(G, tmp, tmp);
					objective = new Conductance(G, tmp, bound);
					trimming = new BoundarySharpness();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (first == "Conductance"
						&& second == "NodeClusterSimilarity"
						&& third == "Dummy") {
					similarity = new NodeClusterSimilarity(G, tmp, tmp);
					objective = new Conductance(G, tmp, bound);
					trimming = new DummyTrimming();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (first == "Conductance" && second == "Dummy"
						&& third == "BoundarySharpness") {
					similarity = new DummySimilarity(G, tmp, tmp);
					objective = new Conductance(G, tmp, bound);
					trimming = new BoundarySharpness();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (first == "Conductance" && second == "Dummy"
						&& third == "Dummy") {
					similarity = new DummySimilarity(G, tmp, tmp);
					objective = new Conductance(G, tmp, bound);
					trimming = new DummyTrimming();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (first == "ModularityM"
						&& second == "NodeClusterSimilarity"
						&& third == "BoundarySharpness") {
					similarity = new NodeClusterSimilarity(G, tmp, tmp);
					objective = new LocalModularityM(G, tmp, bound);
					trimming = new BoundarySharpness();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (first == "ModularityM"
						&& second == "NodeClusterSimilarity"
						&& third == "Dummy") {
					similarity = new NodeClusterSimilarity(G, tmp, tmp);
					objective = new LocalModularityM(G, tmp, bound);
					trimming = new DummyTrimming();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (first == "ModularityM" && second == "Dummy"
						&& third == "BoundarySharpness") {
					similarity = new DummySimilarity(G, tmp, tmp);
					objective = new LocalModularityM(G, tmp, bound);
					trimming = new BoundarySharpness();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (first == "ModularityM" && second == "Dummy"
						&& third == "Dummy") {
					similarity = new DummySimilarity(G, tmp, tmp);
					objective = new LocalModularityM(G, tmp, bound);
					trimming = new DummyTrimming();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (first == "ModularityL"
						&& second == "NodeClusterSimilarity"
						&& third == "BoundarySharpness") {
					similarity = new NodeClusterSimilarity(G, tmp, tmp);
					objective = new LocalModularityL(G, tmp, bound);
					trimming = new BoundarySharpness();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (first == "ModularityL"
						&& second == "NodeClusterSimilarity"
						&& third == "Dummy") {
					similarity = new NodeClusterSimilarity(G, tmp, tmp);
					objective = new LocalModularityL(G, tmp, bound);
					trimming = new DummyTrimming();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (first == "ModularityL" && second == "Dummy"
						&& third == "BoundarySharpness") {
					similarity = new DummySimilarity(G, tmp, tmp);
					objective = new LocalModularityL(G, tmp, bound);
					trimming = new BoundarySharpness();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (first == "ModularityL" && second == "Dummy"
						&& third == "Dummy") {
					similarity = new DummySimilarity(G, tmp, tmp);
					objective = new LocalModularityL(G, tmp, bound);
					trimming = new DummyTrimming();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else {
					std::cout << "[ERROR] invalid arguments " << std::endl;
					exit(1);
				}

			} else if (Aux::StringTools::split(detectorArg, ':').size() == 3) {
				std::string first = Aux::StringTools::split(detectorArg, ':')[1];
				std::string second =
						Aux::StringTools::split(detectorArg, ':')[2];
				if (first == "Conductance") {
					objective = new Conductance(G, tmp, bound);
					if (second == "NodeClusterSimilarity") {
						similarity = new NodeClusterSimilarity(G, tmp, tmp);
						trimming = new DummyTrimming();
						algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
					} else if (second == "BoundarySharpness") {
						similarity = new DummySimilarity(G, tmp, tmp);
						trimming = new BoundarySharpness();
						algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
					} else if (second == "Dummy") {
						similarity = new DummySimilarity(G, tmp, tmp);
						trimming = new DummyTrimming();
						algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
					} else {
						std::cout << "[ERROR] invalid arguments " << std::endl;
						exit(1);
					}
				} else if (first == "ModularityM") {
					objective = new LocalModularityM(G, tmp, bound);
					if (second == "NodeClusterSimilarity") {
						similarity = new NodeClusterSimilarity(G, tmp, tmp);
						trimming = new DummyTrimming();
						algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
					} else if (second == "BoundarySharpness") {
						similarity = new DummySimilarity(G, tmp, tmp);
						trimming = new BoundarySharpness();
						algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
					} else if (second == "Dummy") {
						similarity = new DummySimilarity(G, tmp, tmp);
						trimming = new DummyTrimming();
						algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
					} else {
						std::cout << "[ERROR] invalid arguments " << std::endl;
						exit(1);
					}
				} else if (first == "ModularityL") {
					objective = new LocalModularityL(G, tmp, bound);
					if (second == "NodeClusterSimilarity") {
						similarity = new NodeClusterSimilarity(G, tmp, tmp);
						trimming = new DummyTrimming();
						algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
					} else if (second == "BoundarySharpness") {
						similarity = new DummySimilarity(G, tmp, tmp);
						trimming = new BoundarySharpness();
						algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);;
					} else if (second == "Dummy") {
						similarity = new DummySimilarity(G, tmp, tmp);
						trimming = new DummyTrimming();
						algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
					} else {
						std::cout << "[ERROR] invalid arguments " << std::endl;
						exit(1);
					}
				} else {
					std::cout << "[ERROR] invalid arguments " << std::endl;
					exit(1);
				}
			} else if (Aux::StringTools::split(detectorArg, ':').size() == 2) {
				if (Aux::StringTools::split(detectorArg, ':')[1]
						== "Conductance") {
					objective = new Conductance(G, tmp, bound);
					similarity = new DummySimilarity(G, tmp, tmp);
					trimming = new DummyTrimming();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (Aux::StringTools::split(detectorArg, ':')[1]
						== "ModularityM") {
					objective = new LocalModularityM(G, tmp, bound);
					similarity = new DummySimilarity(G, tmp, tmp);
					trimming = new DummyTrimming();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else if (Aux::StringTools::split(detectorArg, ':')[1]
						== "ModularityL") {
					objective = new LocalModularityL(G, tmp, bound);
					similarity = new DummySimilarity(G, tmp, tmp);
					trimming = new DummyTrimming();
					algo = new GreedyCommunityExpansion(G, *similarity, *objective, *trimming);
				} else {
					std::cout << "[ERROR] unknown objective: " << Aux::StringTools::split(detectorArg, ':')[0] << std::endl;
					exit(1);
				}
			} else {
				std::cout << "[ERROR] invalid argument " << std::endl;
				exit(1);
			}
		}
	} else {
		std::cout
				<< "[ERROR] option --detector=<NAME>:<PARAMS> must be supplied"
				<< std::endl;
		exit(1);
	}
	std::cout << "[BEGIN]" << std::endl;


	assert (algo != NULL);


	// RUN
	Aux::Timer running1;
	Aux::Timer running2;

	std::unordered_map<count, int64_t> timeMap;
	std::unordered_map<count, std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>>> results;
	int64_t runtime;

	running1.start();
	for (count i = 0; i < runs; i++) {
		running2.start();
		std::unordered_set<node> seeds = seedGen->getSeeds(nSeeds);
		std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> result = algo->run(
				seeds);
		running2.stop();
		results.insert({i, result});
		timeMap.insert({i, running2.elapsedMilliseconds()});
	}
	running1.stop();
	runtime = running1.elapsedMilliseconds();
	std::cout << "[DONE]" << std::endl;


	// EVALUATION
	if (options[SUMMARY]) {
		std::ofstream summary(options[SUMMARY].arg);
		summary << "Node ID" << ";" << "Conductance" << ";" << "Local Modularity" << ";"
				<< "Community Size" << ";" << "Runtime" << std::endl;
		for (auto u : results) {
			for (auto v : u.second) {
				if (measures.size() == 0) {
					summary << v.first << ";" << v.second.second << std::endl;
				} else if (measures.size() == 1) {
					summary << v.first << ";"
							<< (measures[0])->getQuality(v.second.first) << ";"
							<< v.second.first.size() << ";"
							<< v.second.second << std::endl;
				} else if (measures.size() == 2) {
					summary << v.first << ";"
							<< (measures[0])->getQuality(v.second.first) << ";"
							<< (measures[1])->getQuality(v.second.first) << ";"
							<< v.second.first.size() << ";"
							<< v.second.second << std::endl;
				}
			}
		}
	}
	std::unordered_set<node> seed;
	if (options[SHOW]) {
		std::string showArg = options[SHOW].arg;
		if (Aux::StringTools::split(showArg, ':').size() == 2) {
			std::string seedNode = Aux::StringTools::split(showArg, ':')[0];
			std::string path = Aux::StringTools::split(showArg, ':')[1];
			node tmp = std::stoi(seedNode);
			seed.insert(tmp);
			std::unordered_set<node> community = algo->run(seed).find(tmp)->second.first;
			Graph sub = Subgraph::fromNodes(G, community);
			METISGraphWriter writer;
			writer.write(sub, path);
		}
	}
	// optionally get ground truth
	if (options[GROUND_TRUTH]) {
		// TODO: @Alex - get LFR ground truth
		// TODO: evaluate ground truth
	} else {
		std::cout << "[INFO]�no ground truth supplied" << std::endl;
	}

	std::cout << "[EXIT] terminated normally" << std::endl;
	return 0;
}



