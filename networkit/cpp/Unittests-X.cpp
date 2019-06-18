/*
 * Unittests-X.cpp
 *
 *  Created on: 2012
 *      Author: Christian Staudt <christian.staudt@kit.edu>,
 *              Henning Meyerhenke <henning.meyerhenke@kit.edu>,
 *              Manuel Penschuck <networkit@manuel.jetzt>
 */


// includes
#include <algorithm>
#include <iostream>
#include <utility>

// GoogleTest
#include <gtest/gtest.h>

// OpenMP
#include <omp.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/ext/optionparser.hpp>

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

enum  optionIndex { UNKNOWN, HELP, LOGLEVEL, THREADS, TESTS, RUNNABLE, DEBUG, BENCHMARKS };
const OptionParser::Descriptor usage[] = {
    {UNKNOWN,	0,	"" ,	"",				OptionParser::Arg::None,	"Options:" },
    {HELP,		0,	"h",	"help",			OptionParser::Arg::None,	"  --help  \t Print usage and exit." },
    {LOGLEVEL,	0,	"",		"loglevel",		OptionParser::Arg::Required,"  --loglevel=<LEVEL>  \t set the log level" },
    {THREADS,	0,	"",		"threads",		OptionParser::Arg::Required,"  --threads=<NUM>  \t set the maximum number of threads" },
    {TESTS,		0,	"t",	"tests",		OptionParser::Arg::None,	"  --tests \t Run unit tests"},
    {RUNNABLE,	0,	"r",	"run",			OptionParser::Arg::None,	"  --run \t Run unit tests which don't use assertions"},
    {DEBUG,		0,	"d",	"debug",		OptionParser::Arg::None,	"  --debug \t Run tests to debug some algorithms"},
    {BENCHMARKS,0,	"b",	"benchmarks",	OptionParser::Arg::None,	"  --benchmarks \t Run benchmarks"},
    {UNKNOWN,	0,	"",		"",				OptionParser::Arg::None,	"" },
    {0,			0,	0,		0,				0,							0}
};


int main(int argc, char **argv) {
	std::cout << "*** NetworKit Unit Tests ***\n";

	// PARSE OPTIONS
	const auto parse_argc = std::max(0, argc - 1);
	const auto parse_argv = argv + (argc>0); // skip program name argv[0] if present
	OptionParser::Stats  stats(usage, parse_argc, parse_argv);
	std::vector<OptionParser::Option> options(stats.options_max), buffer(stats.buffer_max);
	OptionParser::Parser parse(usage, parse_argc, parse_argv, options.data(), buffer.data());

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
	omp_set_nested(1); // enable nested parallelism

	if (options[THREADS]) {
		// set number of threads
		int nThreads = std::atoi(options[THREADS].arg);
		Aux::setNumberOfThreads(nThreads);
	}

	::testing::InitGoogleTest(&argc, argv);
	auto set_filter = [&] (const char* filter) {
	    std::cout << "Set --gtest_filter='" << filter << "'\n";
        ::testing::GTEST_FLAG(filter) = filter;
	};

	if (options[TESTS]) {
        set_filter("*Test.test*");
	} else if (options[DEBUG]) {
        set_filter("*Test.debug*");
	} else if (options[BENCHMARKS]) {
        set_filter("*Benchmark*");
	}	else if (options[RUNNABLE]) {
        set_filter("*Test.run*");
	}

	INFO("=== starting unit tests ===");

	return RUN_ALL_TESTS();
}
