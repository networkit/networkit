/*
 * Unittests-X.cpp
 *
 *  Created on: 2012
 *      Author: Christian Staudt <christian.staudt@kit.edu>,
 *              Henning Meyerhenke <henning.meyerhenke@kit.edu>,
 *              Manuel Penschuck <networkit@manuel.jetzt>
 */

#include <algorithm>
#include <iostream>
#include <omp.h>
#include <utility>

#include <gtest/gtest.h>

#include <tlx/cmdline_parser.hpp>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallelism.hpp>

struct Options {
    std::string loglevel;
    unsigned numThreads{0};

    bool modeTests{false};
    bool modeDebug{false};
    bool modeBenchmarks{false};
    bool modeRunnable{false};

    bool parse(int argc, char* argv[]) {
        tlx::CmdlineParser parser;

        parser.add_bool('t', "tests",      modeTests,      "Run unit tests");
        parser.add_bool('r', "run",        modeRunnable,   "Run unit tests which don't use assertions");
        parser.add_bool('d', "debug",      modeDebug,      "Run tests to debug some algorithms");
        parser.add_bool('b', "benchmarks", modeBenchmarks, "Run benchmarks");
        
        parser.add_unsigned("threads",     numThreads,     "set the maximum number of threads; 0 (=default) uses OMP default");
        parser.add_string("loglevel",      loglevel,       "set the log level (TRACE|DEBUG|INFO|WARN|ERROR|FATAL)");

        if (!parser.process(argc, argv, std::cerr))
            return false;

        if (modeTests + modeDebug + modeBenchmarks + modeRunnable > 1) {
            std::cerr << "Select at most one of -t, -r, -d, -b\n";
            return false;
        }

        return true;
    }
};

int main(int argc, char *argv[]) {
    std::cout << "*** NetworKit Unit Tests ***\n";

    ::testing::InitGoogleTest(&argc, argv);
    Options options;
    if (!options.parse(argc, argv)) {
        return -1;
    }

    // Configure logging
#ifndef NOLOGGING
    if (!options.loglevel.empty()) {
        Aux::Log::setLogLevel(options.loglevel);
        if(Aux::Log::getLogLevel() == "INFO"){
            Aux::Log::Settings::setPrintLocation(false);
        }else{
            Aux::Log::Settings::setPrintLocation(true);
        }
    } else {
        Aux::Log::setLogLevel("ERROR");	// with default level
        Aux::Log::Settings::setPrintLocation(true);
    }
    std::cout << "Loglevel: " << Aux::Log::getLogLevel() << "\n";
#endif

    // Configure parallelism
    {
        omp_set_nested(1); // enable nested parallelism

        if (options.numThreads) {
            Aux::setNumberOfThreads(options.numThreads);
        }
        std::cout << "Max. number of threads: " << Aux::getMaxNumberOfThreads() << "\n";
    }

    // Configure test filter
    {
        auto setFilter = [&](const char *filter) {
            ::testing::GTEST_FLAG(filter) = filter;
        };

        if (options.modeTests) {
            setFilter("*Test.test*");
        } else if (options.modeDebug) {
            setFilter("*Test.debug*");
        } else if (options.modeBenchmarks) {
            setFilter("*Benchmark*");
        } else if (options.modeRunnable) {
            setFilter("*Test.run*");
        }
    }

    INFO("=== starting unit tests ===");

    return RUN_ALL_TESTS();
}
