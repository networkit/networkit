/*
 * IndependentSetTest.cpp
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/independentset/Luby.hpp>

namespace NetworKit {

class IndependentSetGTest : public testing::Test {};

TEST_F(IndependentSetGTest, testLuby) {
    Aux::Random::setSeed(42, true);
    for (count n = 50; n < 500; n += 50) {
        ErdosRenyiGenerator generator(n, 0.01);
        Graph G = generator.generate();

        Luby luby;
        std::vector<bool> I = luby.run(G);

        EXPECT_TRUE(luby.isIndependentSet(I, G));
    }
}

TEST_F(IndependentSetGTest, debugLuby) {
    Aux::Random::setSeed(42, true);
    count n = 500;
    ErdosRenyiGenerator generator(n, 0.001);
    Graph G = generator.generate();

    Luby luby;
    std::vector<bool> I = luby.run(G);

    EXPECT_TRUE(luby.isIndependentSet(I, G)) << "result must be an independent set";

    count size = 0;
    for (bool x : I) {
        if (x) {
            size += 1;
        }
    }
    INFO("independent set size: ", size, "/", n);
}

TEST_F(IndependentSetGTest, debugLubyWithSelfLoops) {
    Aux::Random::setSeed(42, true);
    count n = 500;
    ErdosRenyiGenerator generator(n, 0.001);
    Graph G = generator.generate();

    G.forNodes([&](node u) { G.addEdge(u, u); });

    Luby luby;
    std::vector<bool> I = luby.run(G);

    EXPECT_TRUE(luby.isIndependentSet(I, G)) << "result must be an independent set";

    count size = 0;
    for (bool x : I) {
        if (x) {
            size += 1;
        }
    }
    INFO("independent set size: ", size, "/", n);
}

} /* namespace NetworKit */
