//============================================================================
// Name        : Unittests-X.cpp
// Author      : Christian Staudt (christian.staudt@kit.edu),
//		 Henning Meyerhenke (henning.meyerhenke@kit.edu)
// Version     :
// Copyright   : � 2012, Christian Staudt, Henning Meyerhenke
// Description : Calling unit tests and benchmarks
//============================================================================

// includes
#include <iostream>
#include <utility>
//#include <cfenv>	// floating point exceptions
#include <stdexcept>

// GoogleTest
#ifndef NOGTEST
#include <gtest/gtest.h>
#endif

// OpenMP
#include <omp.h>

// necessary for some reasons?
#include "Globals.h"
#include "ext/optionparser.h"
#include "auxiliary/Log.h"
#include "graph/Graph.h"
#include "auxiliary/Parallelism.h"


using namespace NetworKit;


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

// TODO: clean up obsolete parameters
enum  optionIndex { UNKNOWN, HELP, LOGLEVEL, THREADS, TESTS, RUNNABLE, DEBUG, BENCHMARKS, FILTER };
const OptionParser::Descriptor usage[] =
{
{UNKNOWN,	0,	"" ,	"",				OptionParser::Arg::None,	"Options:" },
{HELP,		0,	"h",	"help",			OptionParser::Arg::None,	"  --help  \t Print usage and exit." },
{LOGLEVEL,	0,	"",		"loglevel",		OptionParser::Arg::Required,"  --loglevel=<LEVEL>  \t set the log level" },
{THREADS,	0,	"",		"threads",		OptionParser::Arg::Required,"  --threads=<NUM>  \t set the maximum number of threads" },
{TESTS,		0,	"t",	"tests",		OptionParser::Arg::None,	"  --tests \t Run unit tests"},
{RUNNABLE,	0,	"r",	"run",			OptionParser::Arg::None,	"  --run \t Run unit tests which don't use assertions"},
{DEBUG,		0,	"d",	"debug",		OptionParser::Arg::None,	"  --debug \t Run tests to debug some algorithms"},
{BENCHMARKS,0,	"b",	"benchmarks",	OptionParser::Arg::None,	"  --benchmarks \t Run benchmarks"},
{FILTER,	0,	"f",	"gtest_filter",	OptionParser::Arg::Required,"  --gtest_filter=<FILTER_PATTERN> \t Run tests that match the filter pattern" },
{UNKNOWN,	0,	"",		"",				OptionParser::Arg::None,	"\nExamples:\n TODO" },
{0,			0,	0,		0,				0,							0}
};


int main(int argc, char **argv) {
	std::cout << "*** NetworKit Unit Tests *** " << std::endl;

	// ENABLE FLOATING POINT EXCEPTIONS (needs GNU extension, apparently only available on Linux)
#ifdef _GNU_SOURCE
	// feenableexcept(FE_ALL_EXCEPT);
#endif

	std::string program_name = argv[0];
	// PARSE OPTIONS
	argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present

	OptionParser::Stats  stats(usage, argc, argv);
	std::vector<OptionParser::Option> options(stats.options_max), buffer(stats.buffer_max);
	OptionParser::Parser parse(usage, argc, argv, options.data(), buffer.data());

	if (parse.error())
	 return 1;

	if (options[HELP]) {
	 OptionParser::printUsage(std::cout, usage);
	 return 0;
	}

	for (OptionParser::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
	 std::cout << "Unknown option: " << opt->name << "\n";

	for (int i = 0; i < parse.nonOptionsCount(); ++i)
	 std::cout << "Non-option #" << i << ": " << parse.nonOption(i) << "\n";




	// CONFIGURE LOGGING
#ifndef NOLOGGING
	if (options[LOGLEVEL]) {
		Aux::Log::setLogLevel(options[LOGLEVEL].arg);
		if(Aux::Log::getLogLevel() == "INFO"){
			Aux::Log::Settings::setPrintLocation(false);
		}else{
			Aux::Log::Settings::setPrintLocation(true);
		}
	} else {
		Aux::Log::setLogLevel("ERROR");	// with default level
		Aux::Log::Settings::setPrintLocation(true);
	}
#endif


	// CONFIGURE PARALLELISM
#ifdef _OPENMP
	omp_set_nested(1); // enable nested parallelism
#endif

	if (options[THREADS]) {
		// set number of threads
		int nThreads = std::atoi(options[THREADS].arg);
		Aux::setNumberOfThreads(nThreads);
	}

	// get program name (currently only for unix)
	auto pos = program_name.find_last_of("/");
	program_name = program_name.substr(pos+1,program_name.length()-1);

#ifndef NOGTEST
	if (options[TESTS]) {
		::testing::GTEST_FLAG(filter) = "*Test.test*";
	} else if (options[DEBUG]) {
		::testing::GTEST_FLAG(filter) = "*Test.debug*";
	} else if (options[BENCHMARKS]) {
		if (program_name != "NetworKit-Tests-O") {
			std::cout << "Hint: Performance tests should be run in optimized mode" << std::endl;
		}
		::testing::GTEST_FLAG(filter) = "*Benchmark*";
	} else if (options[FILTER]) {
		::testing::GTEST_FLAG(filter) = options[FILTER].arg;
	}	else if (options[RUNNABLE]) {
		::testing::GTEST_FLAG(filter) = "*Test.run*";
	}
	::testing::GTEST_FLAG(shuffle) = 1;
	::testing::InitGoogleTest(&argc, argv);
	INFO("=== starting unit tests ===");
	return RUN_ALL_TESTS();
#else
	   throw std::runtime_error("unit tests are excluded from build by the NOGTEST preprocessor directive");
#endif
}
