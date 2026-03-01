/*  EpidemicSimulationSEIRTest.cpp
 *
 *  Created on: 22.02.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <networkit/simulation/EpidemicSimulationSEIR.hpp>
#include <networkit/graph/Graph.hpp>
#include <gtest/gtest.h>

namespace NetworKit {

TEST(TestEpidemicSimulationSEIR, testConstructor){
    Graph graph(1);
    constexpr count maxTimestamps{0};
    constexpr double transmissionProb{0.0};
    constexpr count exposureTime{0};
    constexpr count infectiousTime{0};
    constexpr node startingNode{0};
    EXPECT_NO_THROW(EpidemicSimulationSEIR simulator(graph,maxTimestamps,transmissionProb,exposureTime,infectiousTime,startingNode));
}

TEST(TestEpidemicSimulationSEIR, testNoMaxTimestampsNoStats){
    Graph graph(2);
    graph.addEdge(0,1);
    constexpr count maxTimestamps{0};
    constexpr double transmissionProb{1.0};
    constexpr count exposureTime{2};
    constexpr count infectiousTime{2};
    constexpr node startingNode{0};
    EpidemicSimulationSEIR simulator(graph,maxTimestamps,transmissionProb,exposureTime,infectiousTime,startingNode);
    simulator.run();
    auto stats = simulator.getData();
    EXPECT_TRUE(stats.empty());
}

TEST(TestEpidemicSimulationSEIR, testZeroTransmissionProbabalityNoSpread){
    Graph graph(10);
    graph.addEdge(0,1);
    constexpr count states{4};
    constexpr count maxTimestamps{10};
    constexpr double transmissionProb{0.0};
    constexpr count exposureTime{2};
    constexpr count infectiousTime{2};
    constexpr node startingNode{0};
    EpidemicSimulationSEIR simulator(graph,maxTimestamps,transmissionProb,exposureTime,infectiousTime,startingNode);
    simulator.run();
    auto stats = simulator.getData();
    EXPECT_EQ(stats.size(), maxTimestamps*states);
    auto popAt = [&](count t, count s) -> count {
        const count base = t * states;
        const count idx = base + (s - 0); // assumes State::S == 0
        return stats.at(idx).at(3);
    };

    for (count time = 0; time < maxTimestamps; ++time) {
        const count Susceptible = popAt(time, 0);
        const count Exposed = popAt(time, 1);
        const count Infectious = popAt(time, 2);
        const count Removed = popAt(time, 3);

        EXPECT_EQ(Exposed, 0);
        EXPECT_EQ(Susceptible, 9);
        if (time - infectiousTime >= 0)
            EXPECT_EQ(Removed, 1);
        else
            EXPECT_EQ(Infectious, 1);
        EXPECT_EQ(Susceptible + Exposed + Infectious + Removed, 10);
    }
}
}

